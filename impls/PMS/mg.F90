!-----------------------------------------------------------------
!     Ravi Samtaney & Mark Adams
!     Copyright 2014
!-----------------------------------------------------------------
subroutine solve(g,Apply1,Apply2,Relax1,Res0,errors,nViters)
  use GridModule
  use domain, only:dom_max_grids,nvcycles,nfcycles,rtol,verbose
  use mpistuff
  use error_data_module
  implicit none
  !
  type(proc_patch),intent(in)::g(0:dom_max_grids-1)
  type(error_data),intent(out)::errors(0:dom_max_grids-1+nvcycles)
  external Apply1, Apply2, Relax1
  double precision,intent(out)::Res0
  integer,intent(out):: nViters
  !     Local Vars
  integer:: iter,coarsest_grid,ii
  double precision:: Res,Reslast
  double precision,dimension(&
       g(0)%ilo:g(0)%ihi, g(0)%jlo:g(0)%jhi,&
       g(0)%klo:g(0)%khi, nvar) :: L2u_f,rhs,ux
  double precision norm
  !

  !call Apply2(L2u_f,ux,g(0)) ! ux==0 
  call formRHS(rhs,g(0))
  Res0 = norm(rhs,g(0),2)
  if(mype==0.and.verbose.gt.0)write(6,'(A,E14.6)') '0) solve: |f|=',Res0
  Reslast = Res0
  ! FMG solve
  iter = 0
  err_lev = 0
  if (nfcycles .ne. 0) then
     iter = iter + 1
     call MGF(ux,rhs,g,0,Apply1,Apply2,Relax1,errors)
     call Apply2(L2u_f,ux,g(0))
     Res = norm(rhs-L2u_f,g(0),2)
     if(mype==0.and.verbose.gt.0)write(6,'(I2,A,E14.6,A,F10.6)') &
          iter,') solve: |f-Au|_2=',Res,', rate=',Res/Reslast
     Reslast=Res
  else
     ! initailize what is done in FMG
     call formRHS(rhs,g(0)) ! form RHS explicity on fine grid for vycycles
     Res = Res0
  end if
  
  ! find coarsest grid for diagnostics
  coarsest_grid = -1
  do ii=0,dom_max_grids-1
     if( is_top(g(ii)) ) then
        coarsest_grid = ii
        exit
     end if
  end do

  ! finish with V-cycle
  nViters = 0
  do while(Res/Res0 > rtol .and. iter < nvcycles)
     !call Apply1(L1u,ux,g(0))
     !L2u_f = rhs - L2u_f + L1u ! defect correction
     !call MGV(ux,L2u_f,g,0,Apply1,Relax1)
     call MGV(ux,rhs,g,0,Apply1,Relax1)
     iter=iter+1 ! 2 for FMG and 1 for V 
     ! residual
     call Apply2(L2u_f,ux,g(0))            
     Res = norm(rhs-L2u_f,g(0),2)
     if(mype==0.and.verbose.gt.0)write(6,'(I2,A,E14.6,A,F10.6)') &
          iter,') solve: |f-Au|_2=',Res,', rate=',Res/Reslast
     
     ! form errors & convergance measure
     errors(err_lev)%resid = Res
     call formExactU(L2u_f,g(0))
     errors(err_lev)%uerror = norm(L2u_f-ux,g(0),error_norm)
     call formGradError(L2u_f,ux,g(0),errors(err_lev)%graduerr,error_norm)
     err_lev = err_lev + 1

     Reslast=Res
     nViters = nViters + 1
  enddo
  return
end subroutine solve
!-----------------------------------------------------------------------
recursive subroutine MGF(ux,rhs,g,lev,Apply1,Apply2,Relax1,errors)
  use GridModule
  use mpistuff,only:mype
  use error_data_module
  use domain, only:dom_max_grids,nfmgvcycles,ncoarsesolveits,verbose,nvcycles
  implicit none
  integer,intent(in) :: lev
  type(proc_patch),intent(in)::g(0:dom_max_grids-1)
  double precision,intent(in),dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: rhs
  double precision,intent(out),dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: ux
  external Apply1, Apply2, Relax1
  type(error_data),intent(out)::errors(0:dom_max_grids-1+nvcycles)
  ! Local Vars - Coarse. This is new data on the stack
  integer :: ii,jj,kk
  double precision :: norm,tt
  ! Coarse data. This is new data on the stack
  double precision,dimension(g(lev+1)%ilo:g(lev+1)%ihi,g(lev+1)%jlo:g(lev+1)%jhi,&
       g(lev+1)%klo:g(lev+1)%khi, nvar):: uxC,FC
  double precision,dimension(g(lev)%ilo:g(lev)%ihi,g(lev)%jlo:g(lev)%jhi,&
       g(lev)%klo:g(lev)%khi, nvar)::tmp

  ux = 0.d0 ! everthing uses initial solution except FMG
  if( is_top(g(lev)) ) then
     ! the real start of FMG - coarse grid solve
     if (lev.gt.0 .and. is_top(g(lev-1))) stop 'MGF: not coarsest grid'
     if(mype==0.and.verbose.gt.4)write(6,'(I2,A,I2,A,E14.6)') '[ 0]MGF: bot solve lev=',lev
     !call MGV(ux,rhs,g,lev,Apply1,Relax1)
     call Relax1(ux,rhs,g(lev),ncoarsesolveits)  ! coar grid solve
     if (verbose.gt.1) then
        tt = norm(ux,g(lev),2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') MGF bot |u_1|=',tt
     end if
  else
     ! go "down" the V, allating data (on stack), forming RHS, zero out u
     call formRHS(FC,g(lev+1)) ! form RHS explicity on fine grid
     call MGF(uxC,FC,g,lev+1,Apply1,Apply2,Relax1,errors) 
     ! prolongate and correct
     call Prolong_2(ux,uxC,g(lev),g(lev+1))
     if (verbose.gt.1) then
        tt = norm(ux,g(lev),2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') MGF |u_1|=',tt
     end if
     ! start of multiple v-cycle
     jj = nfmgvcycles !; if(lev==0) jj=10
     do ii=1,jj
        call MGV(ux,rhs,g,lev,Apply1,Relax1)
     enddo
     if (verbose.gt.1.and.mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev
  end if

  ! form errors & convergance measure, recursive so goes from coarse to fine (good)
  call formExactU(tmp,g(lev))
  errors(err_lev)%uerror = norm(tmp-ux,g(lev),error_norm)
  call formGradError(tmp,ux,g(lev),errors(err_lev)%graduerr,error_norm)
  call Apply2(tmp,ux,g(lev))
  errors(err_lev)%resid = norm(tmp-rhs,g(lev),2)
 
  if (verbose.gt.0) then
     if(mype==0)write(6,'(A,I2,A,E14.6,A,E14.6,A,I8,A,I4)') &
          '     lev=',lev,') MGF |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror,&
          ', nx=',g(lev)%imax*g(lev)%nprocx,', npe x=',g(lev)%nprocx
  end if
  
  err_lev = err_lev + 1
  return
end subroutine MGF
!-----------------------------------------------------------------------
recursive subroutine MGV(ux,rhs,g,lev,Apply,Relax)
  use GridModule
  use mpistuff, only:mype
  use domain, only:dom_max_grids,verbose,ncoarsesolveits,ncycles,nsmoothsdown,nsmoothsup,mg_min_size
  !
  implicit none
  integer,intent(in):: lev
  type(proc_patch),intent(in)::g(0:dom_max_grids-1)
  external Apply,Relax
  !
  double precision,dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi, &
       g(lev)%klo:g(lev)%khi, nvar) :: ux, rhs
  double precision norm,tt
  integer:: ii,jj
  !     Coarse 
  double precision,dimension(&
       g(lev+1)%ilo:g(lev+1)%ihi, g(lev+1)%jlo:g(lev+1)%jhi, &
       g(lev+1)%klo:g(lev+1)%khi,nvar) :: uxC, tmpC, FC, ResC
  !     This level
  double precision,dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: Res
  if (verbose.gt.1) then  
     tt = norm(ux,g(lev),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: u_0=',tt
     tt = norm(rhs,g(lev),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: f_0=',tt
  end if
  if( is_top(g(lev)) ) then
     call Relax(ux,rhs,g(lev),ncoarsesolveits)
  else
     !     pre smoothing
     call Relax(ux,rhs,g(lev),nsmoothsdown)
     if (verbose.gt.2) then
        tt = norm(ux,g(lev),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: u_1=',tt
     end if
     !     restrict residual
     call Apply(Res,ux,g(lev))
     Res = rhs - Res
     if (verbose.gt.2) then
        tt = norm(Res,g(lev),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: res^f=',tt
     end if
     call Restrict(FC,Res,g(lev+1),g(lev),1)
     if (verbose.gt.2) then
        tt = norm(FC,g(lev+1),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: res^c_0=',tt
     end if
     !     restrict solution
     call Restrict(uxC,ux,g(lev+1),g(lev),1)
     !     rhs = residual + Ac(Uc)
     call Apply(tmpC,uxC,g(lev+1))
     FC = FC + tmpC
     if (verbose.gt.2) then
        tt = norm(FC,g(lev+1),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: res^c_1=',tt
     end if
     !     Temporarily store uxC into tmpC
     tmpC = uxC
     !     start of v(1) or w(2) cycle 
     jj = ncycles ! ; if(lev==0) jj=40 ! two level method, sort of
     do ii=1,jj
        call MGV(uxC,FC,g,lev+1,Apply,Relax)
     enddo
     !     subtract off old Uc
     uxC = uxC - tmpC
     ! prolongate and correct
     call Prolong_2(ux,uxC,g(lev),g(lev+1))
     if (verbose.gt.2) then
        tt = norm(ux,g(lev),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: u_2=',tt
     end if
     ! end of v cycle, post smoothing
     call Relax(ux,rhs,g(lev),nsmoothsup)
     if (verbose.gt.2) then
        tt = norm(ux,g(lev),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: u_3=',tt
     end if
  end if
  return
end subroutine MGV
!-----------------------------------------------------------------
subroutine destroy_grids_private(gg,max_sz)
  use GridModule
  use mpistuff
  use domain, only:mg_min_size
  implicit none
  integer,intent(in):: max_sz
  type(proc_patch):: gg(0:max_sz-1)
  ! 
  integer :: n
  do n=1,max_sz-1
     !       take care of reductions       
     if (is_top(gg(n-1))) then
        exit
     else if (gg(n-1)%imax > mg_min_size .and. gg(n-1)%jmax > mg_min_size &
#ifdef TWO_D
          ) then
#else
        .and. gg(n-1)%kmax > mg_min_size ) then
#endif
        ! normal reduction, no split yet -- nothing to destroy
     else 
        ! reduce number of procs on grid by 2
        if (gg(n)%nprocx .ne. gg(n-1)%nprocx ) then
           call MPI_COMM_free(gg(n)%comm,ierr)
           call MPI_COMM_free(gg(n)%loc_comm,ierr) 
           call MPI_comm_free(gg(n)%comm3d, ierr)
        end if
     end if
     !       
  enddo
  return
  end subroutine destroy_grids_private
!-----------------------------------------------------------------
subroutine new_grids_private(gg,NProcAxis,iProcAxis,nx,ny,nz,&
     nxlocal,nylocal,nzlocal,comm3d)
  use GridModule
  use mpistuff
  use domain
  implicit none
  integer,intent(in):: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
  integer :: NProcAxis(3),iProcAxis(3)
  type(proc_patch),intent(out):: gg(0:dom_max_grids-1)
  integer :: n,ndims,ii,jj,kk,proc_id,rank
  logical :: periods(3)
  ndims = 3

  gg(0)%dxg=(xr-xl)/nx
  gg(0)%dyg=(yr-yl)/ny
  gg(0)%dzg=(zr-zl)/nz

  gg(0)%comm3d = comm3d
  gg(0)%comm = MPI_COMM_WORLD
  gg(0)%loc_comm = mpi_comm_null 
  !     left, ...
  call MPI_Cart_Shift(comm3D,0,1,gg(0)%left,gg(0)%right,ierr)
  call MPI_Cart_Shift(comm3D,1,1,gg(0)%bottom,gg(0)%top,ierr)
  call MPI_cart_Shift(comm3D,2,1,gg(0)%behind,gg(0)%forward,ierr)
  if(mype==0.and.verbose.gt.4) write(6,*)'[',mype,'] l,r,t,b=',&
       gg(0)%left,gg(0)%right,gg(0)%top,&
       gg(0)%bottom,gg(0)%behind,gg(0)%forward
  !       iprocx ...	
  gg(0)%iprocx = iProcAxis(1)
  gg(0)%iprocy = iProcAxis(2)
  gg(0)%iprocz = iProcAxis(3)
  
  gg(0)%nprocx = NProcAxis(1) 
  gg(0)%nprocy = NProcAxis(2) 
  gg(0)%nprocz = NProcAxis(3) 
  ! derived quantities
  gg(0)%imax=nxlocal
  gg(0)%jmax=nylocal
  gg(0)%ilo=-ng+1
  gg(0)%ihi=nxlocal+ng
  gg(0)%jlo=-ng+1
  gg(0)%jhi=nylocal+ng
#ifdef TWO_D
  gg(0)%kmax=1
  gg(0)%klo=1
  gg(0)%khi=gg(0)%kmax
#else 
  gg(0)%kmax=nzlocal
  gg(0)%klo=-ng+1
  gg(0)%khi=nzlocal+ng
#endif
  ! cache global indices -- could use these for ilo,...
  gg(0)%iglobalx = getIglobalx(gg(0)) ! one based
  gg(0)%iglobaly = getIglobaly(gg(0))
  gg(0)%iglobalz = getIglobalz(gg(0))
  !
  periods(1) = .false.
  periods(2) = .false.
  periods(3) = .false.
#ifdef XPERIODIC
  periods(1) = .true.
#endif
#ifdef YPERIODIC
  periods(2) = .true.
#endif  
#ifdef ZPERIODIC
  periods(3) = .true.
#endif  
  ! start at 1 as 0 was just done
  do n=1,dom_max_grids-1
     gg(n)%dxg = gg(n-1)%dxg*2.d0
     gg(n)%dyg = gg(n-1)%dyg*2.d0
#ifdef TWO_D
     gg(n)%kmax = 1
     gg(n)%dzg  = gg(n-1)%dzg
#else 
     gg(n)%dzg = gg(n-1)%dzg*2.d0
#endif
     ! take care of reductions
     if( is_top(gg(n-1)) ) then
        ! all done - clear rest of grids
        do ii=n,dom_max_grids-1
           gg(ii)%nprocx = 0 
           gg(ii)%nprocy = 0 
           gg(ii)%nprocz = 0 
           gg(ii)%imax = 0 
           gg(ii)%jmax = 0 
           gg(ii)%kmax = 0 
           ! malloc of next grid in MGV, etc. put some dummy
           gg(ii)%ilo=0
           gg(ii)%ihi=1
           gg(ii)%jlo=0
           gg(ii)%jhi=1
           gg(ii)%klo=0
           gg(ii)%khi=1
        end do
        exit
     else if (gg(n-1)%imax>mg_min_size .and. gg(n-1)%jmax>mg_min_size &
#ifndef TWO_D
          .and. gg(n-1)%kmax>mg_min_size &
#endif
          ) then
        !     normal reduction, no split yet	      
        gg(n)%imax=gg(n-1)%imax/2
        gg(n)%jmax=gg(n-1)%jmax/2
#ifndef TWO_D
        gg(n)%kmax=gg(n-1)%kmax/2
#else
        gg(n)%kmax=gg(n-1)%kmax ! 1
#endif
        gg(n)%comm3d = comm3d
        gg(n)%comm = MPI_COMM_WORLD
        gg(n)%iprocx = gg(n-1)%iprocx
        gg(n)%iprocy = gg(n-1)%iprocy
        gg(n)%iprocz = gg(n-1)%iprocz
        gg(n)%nprocx = gg(n-1)%nprocx
        gg(n)%nprocy = gg(n-1)%nprocy
        gg(n)%nprocz = gg(n-1)%nprocz
        
        gg(n)%left = gg(n-1)%left
        gg(n)%right = gg(n-1)%right
        gg(n)%top = gg(n-1)%top
        gg(n)%bottom = gg(n-1)%bottom
        gg(n)%forward = gg(n-1)%forward
        gg(n)%behind = gg(n-1)%behind
        gg(n)%loc_comm = mpi_comm_null
     else  ! gg(n-1)%imax == min_psize ...
        !     reduce number of procs on grid by 2
        if(gg(n-1)%nprocz==1 .or. gg(n-1)%nprocy==1 &
#ifndef TWO_D
             .or. gg(n-1)%nprocx==1 ) then ! one proc reduction, copy comm stuff
#else
           ) then
#endif
           gg(n)%nprocx = gg(n-1)%nprocx
           gg(n)%nprocy = gg(n-1)%nprocy
           gg(n)%iprocx = gg(n-1)%iprocx
           gg(n)%iprocy = gg(n-1)%iprocy
           gg(n)%nprocz = gg(n-1)%nprocz
           gg(n)%iprocz = gg(n-1)%iprocz
           
           gg(n)%imax=gg(n-1)%imax/2
           gg(n)%jmax=gg(n-1)%jmax/2
#ifndef TWO_D
           gg(n)%kmax=gg(n-1)%kmax/2
#else
           gg(n)%kmax=gg(n-1)%kmax ! 1
#endif
           ! use old comms, really just comm_self
           gg(n)%comm3d = gg(n-1)%comm3d
           gg(n)%comm = gg(n-1)%comm
           gg(n)%loc_comm = gg(n-1)%loc_comm

           gg(n)%left = gg(n-1)%left
           gg(n)%right = gg(n-1)%right
           gg(n)%top = gg(n-1)%top
           gg(n)%bottom = gg(n-1)%bottom
           gg(n)%forward = gg(n-1)%forward
           gg(n)%behind = gg(n-1)%behind
           if (mype==0.and.verbose.gt.4)then
              write(0,*) '[',mype, '] one pe reduction'
              write(0,*) '[',mype, '] X:',gg(n)%iprocx
              write(0,*) '[',mype, '] Y:',gg(n)%iprocy
              write(0,*) '[',mype, '] Z:',gg(n)%iprocz
              write(0,*) '[',mype, '] nx', gg(n)%nprocx
              write(0,*) '[',mype, '] ny', gg(n)%nprocy
              write(0,*) '[',mype, '] nz', gg(n)%nprocz
           endif
        else ! normal split
           gg(n)%nprocx = gg(n-1)%nprocx/2
           gg(n)%nprocy = gg(n-1)%nprocy/2
           gg(n)%iprocx = (gg(n-1)%iprocx-1)/2 + 1
           gg(n)%iprocy = (gg(n-1)%iprocy-1)/2 + 1
#ifdef TWO_D
           gg(n)%nprocz = gg(n-1)%nprocz
           gg(n)%iprocz = gg(n-1)%iprocz
#else
           gg(n)%nprocz = gg(n-1)%nprocz/2
           gg(n)%iprocz = (gg(n-1)%iprocz-1)/2 + 1
#endif
           ! zero based proc ID
           proc_id = (gg(n)%iprocx-1)*gg(n)%nprocy*gg(n)%nprocz &
                + (gg(n)%iprocy-1)*gg(n)%nprocz + (gg(n)%iprocz-1)
           ii = mod(gg(n-1)%iprocx-1,2)
           jj = mod(gg(n-1)%iprocy-1,2)
           kk = mod(gg(n-1)%iprocz-1,2)
#ifdef TWO_D
           ii = jj + 2*ii           ! local zero based ID
#else 
           ii = kk + 2*jj + 4*ii    ! local zero based ID
#endif
           call MPI_COMM_SPLIT(gg(n-1)%comm,ii,proc_id,gg(n)%comm,ierr)
           call MPI_COMM_SPLIT(gg(n-1)%comm,proc_id,ii,gg(n)%loc_comm,ierr)
           
           NProcAxis(1) = gg(n)%nprocx
           NProcAxis(2) = gg(n)%nprocy
           NProcAxis(3) = gg(n)%nprocz
           
           call MPI_cart_create( gg(n)%comm, ndims, NProcAxis, &
                periods, .false., gg(n)%comm3D, ierr)
           
           ! debug              
           call MPI_comm_rank(gg(n)%comm3D,ii,ierr)
           call MPI_Cart_Coords(gg(n)%comm3D,ii,ndims,NProcAxis,ierr)
           if(gg(n)%iprocx.ne.NProcAxis(1)+1) stop '%iprocx'
           if(gg(n)%iprocy.ne.NProcAxis(2)+1) stop '%iprocy'
           if(gg(n)%iprocz.ne.NProcAxis(3)+1) stop '%iprocz'
           if (mype==npe/2.and.verbose.gt.4)then
              write(0,*) '[',mype, '] cart rank',ii
              write(0,*) '[',mype, '] X:',gg(n)%iprocx,NProcAxis(1)+1
              write(0,*) '[',mype, '] Y:',gg(n)%iprocy,NProcAxis(2)+1
              write(0,*) '[',mype, '] Z:',gg(n)%iprocz,NProcAxis(3)+1
              ii = mod(gg(n-1)%iprocx-1,2)
              jj = mod(gg(n-1)%iprocy-1,2)
              write(0,*) '[',mype,'] ii=',ii,'jj=',jj,'proc_id=',proc_id
              write(0,*) '[',mype, '] nx', gg(n)%nprocx
              write(0,*) '[',mype, '] ny', gg(n)%nprocy
              write(0,*) '[',mype, '] nz', gg(n)%nprocz
           endif
           
           call MPI_Cart_Shift(gg(n)%comm3D,0,1,gg(n)%left,gg(n)%right,ierr)! Neighbors
           call MPI_Cart_Shift(gg(n)%comm3D,1,1,gg(n)%bottom,gg(n)%top,ierr)
           call MPI_cart_Shift(gg(n)%comm3D,2,1,gg(n)%behind,gg(n)%forward,ierr)
           !     logical size of grid stays the same
           gg(n)%imax=gg(n-1)%imax
           gg(n)%jmax=gg(n-1)%jmax
#ifndef TWO_D
           gg(n)%kmax=gg(n-1)%kmax
#endif
        endif
     end if
     if(gg(n)%nprocz==0 .or. gg(n)%nprocy==0 .or. &
          gg(n)%nprocx==0 ) then
        gg(n)%nprocz=0; gg(n)%nprocy=0; gg(n)%nprocx=0
        if(mype==0) write(6,*)'[',mype,'] domain is too thin',gg(n)%nprocx==0,gg(n)%nprocy==0,gg(n)%nprocz==0
        stop 'domain is too thin' ! domain is too thin
     endif
     
     ! add ghosts
     gg(n)%ilo=-ng+1
     gg(n)%ihi= gg(n)%imax+ng
     gg(n)%jlo=-ng+1
     gg(n)%jhi=gg(n)%jmax+ng
#ifdef TWO_D
     gg(n)%klo=1
     gg(n)%khi=1
#else 
     gg(n)%klo=-ng+1
     gg(n)%khi=gg(n)%kmax+ng
#endif
     ! cache global indices -- could use these for ilo,...
     gg(n)%iglobalx = getIglobalx(gg(n)) ! one based, could compute on fly
     gg(n)%iglobaly = getIglobaly(gg(n))
     gg(n)%iglobalz = getIglobalz(gg(n))
     
     ! print topology
     if(mype==0.and.verbose.gt.3) then
        write(6,*) '[',mype,'] level ',n, ', nxl=',gg(n)%imax, ', nxl=',gg(n)%imax
        if (mype==-1) then
           write(6,*) '[',mype,'] ipx=',gg(n)%iprocx, ',ipy=',gg(n)%iprocy,',ipz= ',gg(n)%iprocz
           write(6,*) '[',mype,'] npx=',gg(n)%nprocx, ',npy=',gg(n)%nprocy,',npz=',gg(n)%nprocz
           write(6,*) '[',mype,'] nxl=',gg(n)%imax,   ',nyl=',gg(n)%jmax,  ',nzl=',gg(n)%kmax
        endif
     end if
  enddo
end subroutine new_grids_private
