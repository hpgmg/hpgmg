!-----------------------------------------------------------------
!     Ravi Samtaney & Mark Adams
!     Copyright 2014
!-----------------------------------------------------------------
subroutine solve(g,Apply1,Apply2,Relax1,Res0,errors,nViters)
  use GridModule
  use mpistuff
  use error_data_module
  implicit none
  !
  type(proc_patch),intent(in)::g(0:max_grids-1)
  type(error_data),intent(out)::errors(0:max_grids-1+nvcycles)
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
  if(mype==0)write(6,'(A9,E14.6)') '[ 0] |f|=',Res0
  Reslast = Res0
  ! FMG solve
  iter = 0
  err_lev = 0
  if (nfcycles .ne. 0) then
     iter = iter + 1
     call MGF(ux,rhs,g,0,Apply1,Apply2,Relax1,errors)
     call Apply2(L2u_f,ux,g(0))
     Res = norm(rhs-L2u_f,g(0),2)
     if(mype==0)write(6,'(A1,I2,A10,E14.6,A7,F10.6)') &
          '[',iter,'] |f-Au|_2=',Res,', rate=',Res/Reslast
     Reslast=Res
  else
     ! initailize what is done in FMG
     ux = 0.d0
     call formRHS(rhs,g(0)) ! form RHS explicity on fine grid for vycycles
     Res = Res0
  end if
  
  ! find coarsest grid for diagnostics
  coarsest_grid = -1
  do ii=0,max_grids-1
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
     if(mype==0)write(6,'(A1,I2,A10,E14.6,A7,F10.6)') &
          '[',iter,'] |f-Au|_2=',Res,', rate=',Res/Reslast
     
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
  implicit none
  integer,intent(in) :: lev
  type(proc_patch),intent(in)::g(0:max_grids-1)
  double precision,intent(in),dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: rhs
  double precision,intent(out),dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: ux
  external Apply1, Apply2, Relax1
  type(error_data),intent(out)::errors(0:max_grids-1+nvcycles)
  ! Local Vars - Coarse. This is new data on the stack
  integer :: ii,jj
  double precision :: norm
  ! Coarse data. This is new data on the stack
  double precision,dimension(g(lev+1)%ilo:g(lev+1)%ihi,g(lev+1)%jlo:g(lev+1)%jhi,&
       g(lev+1)%klo:g(lev+1)%khi, nvar):: uxC,FC
  double precision,dimension(g(lev)%ilo:g(lev)%ihi,g(lev)%jlo:g(lev)%jhi,&
       g(lev)%klo:g(lev)%khi, nvar)::tmp

  ux = 0.d0 ! everthing uses initial solution except FMG
  if( is_top(g(lev)) ) then
     ! the real start of FMG - coarse grid solve
     if (lev.gt.0 .and. is_top(g(lev-1))) stop 'MGF: not coarsest grid'
     ! defect correction
     !call Apply2(Res2,ux,g(lev)) 
     !call Apply1(Res,ux,g(lev)) 
     !Res = rhs - Res2 + Res
     call MGV(ux,rhs,g,lev,Apply1,Relax1)
     !call Relax1(ux,rhs,g(lev),16)  ! coar grid solves should be parameter
  else
     ! go "down" the V, allating data (on stack), forming RHS, zero out u
     call formRHS(FC,g(lev+1)) ! form RHS explicity on fine grid
     !uxC = 0.d0                
     call MGF(uxC,FC,g,lev+1,Apply1,Apply2,Relax1,errors) 
     ! prolongate and correct
     call Prolong_2(ux,uxC,g(lev),g(lev+1))
     ! start of multiple v-cycle
     jj = nfmgvcycles !; if(lev==0) jj=10
     !if(mype==0)print *,'FMG V-cycles = ',jj
     do ii=1,jj
        ! defect correction
        ! call Apply2(Res2,ux,g(lev)) 
        ! call Apply1(Res,ux,g(lev))
        ! Res = rhs - Res2 + Res
        !     solve
        call MGV(ux,rhs,g,lev,Apply1,Relax1)
     enddo
  endif

  ! form errors & convergance measure, recursive so goes from coarse to fine (good)
  call formExactU(tmp,g(lev))
  errors(err_lev)%uerror = norm(tmp-ux,g(lev),error_norm)
  call formGradError(tmp,ux,g(lev),errors(err_lev)%graduerr,error_norm)
  call Apply2(tmp,ux,g(lev))
  errors(err_lev)%resid = norm(tmp-rhs,g(lev),2)
  
  err_lev = err_lev + 1

  return
end subroutine MGF
!-----------------------------------------------------------------------
recursive subroutine MGV(ux,rhs,g,lev,Apply,Relax)
  use GridModule
  use mpistuff, only:mype
  !
  implicit none
  integer,intent(in):: lev
  type(proc_patch),intent(in)::g(0:max_grids-1)
  external Apply,Relax
  !
  double precision,dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi, &
       g(lev)%klo:g(lev)%khi, nvar) :: ux, rhs
  double precision norm,tt
  integer:: ii,jj
  integer,parameter::nrelax=1 ! parameter !!!
  !     Coarse 
  double precision,dimension(&
       g(lev+1)%ilo:g(lev+1)%ihi, g(lev+1)%jlo:g(lev+1)%jhi, &
       g(lev+1)%klo:g(lev+1)%khi,nvar) :: uxC, tmpC, FC, ResC
  !     This level
  double precision,dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: Res
  !     
  !     tt = norm2(ux,g(lev)); if(my_pe==0)print*,lev,') V0: ux=',tt
  if( is_top(g(lev)) ) then
     call Relax(ux,rhs,g(lev),16)
     !      tt = norm2(ux,g(lev)); if(my_pe==0)print*,lev,') VT: ux=',tt
  else
     !     pre smoothing
     call Relax(ux,rhs,g(lev),nrelax)
     !     restrict residual
     call Apply(Res,ux,g(lev))
     Res = rhs - Res
     call Restrict(FC,Res,g(lev),g(lev+1),1)
     !     restrict solution
     call Restrict(uxC,ux,g(lev),g(lev+1),1)
     !     rhs = residual + Ac(Uc)
     call Apply(tmpC,uxC,g(lev+1))
     FC = FC + tmpC
     !     Temporarily store uxC into tmpC
     tmpC = uxC
     !     start of v(1) or w(2) cycle 
     jj = ncycles ! ; if(lev==0) jj=40 ! two level method, sort of
     do ii=1,jj
        call MGV(uxC,FC,g,lev+1,Apply,Relax)        
        if(lev==-1)then
           call Apply(ResC,uxC,g(lev+1))            
           tt = norm(FC-ResC,g(lev+1),3)
           if(mype==0)write(6,*) '       [',ii,'] RESIDUAL=',tt
        endif
     enddo
     !     subtract off old Uc
     uxC = uxC - tmpC
     ! prolongate and correct
     call Prolong_2(ux,uxC,g(lev),g(lev+1))
     ! end of v cycle, post smoothing
     call Relax(ux,rhs,g(lev),nrelax)
  endif  
  return
end subroutine MGV
!-----------------------------------------------------------------
subroutine destroy_grids_private(gg,max_sz)
  use GridModule
  use mpistuff
  implicit none
  integer,intent(in):: max_sz
  type(proc_patch):: gg(0:max_sz-1)
  ! 
  integer :: n
  do n=1,max_sz-1
     !       take care of reductions       
     if (is_top(gg(n-1))) then
        exit
     else if (gg(n-1)%imax > min_psize .and. gg(n-1)%jmax > min_psize &
#ifdef TWO_D
          ) then
#else
        .and. gg(n-1)%kmax > min_psize ) then
#endif
        ! normal reduction, no split yet -- nothing to destroy
     else 
        !       reduce number of procs on grid by 2          
        call MPI_COMM_free(gg(n)%comm,ierr)
        call MPI_COMM_free(gg(n)%loc_comm,ierr) 
        call MPI_comm_free(gg(n)%comm3d, ierr)
     endif
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
  type(proc_patch),intent(out):: gg(0:max_grids-1)
  integer :: n,ndims,ii,jj,kk,proc_id,rank
  logical :: periods(3)
  ndims = 3
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
  if(mype==-1) write(6,*)'[',mype,'] l,r,t,b=',gg(0)%left,gg(0)%right,gg(0)%top,&
       gg(0)%bottom,gg(0)%behind,gg(0)%forward
  !       iprocx ...	
  gg(0)%iprocx = iProcAxis(1)
  gg(0)%iprocy = iProcAxis(2)
  gg(0)%iprocz = iProcAxis(3)
  !	write (0,*) '[',mype,'] g iprocx=',iprocx,iprocy,iprocz
  gg(0)%xprocs = NProcAxis(1) 
  gg(0)%yprocs = NProcAxis(2) 
  gg(0)%zprocs = NProcAxis(3) 
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
  do n=1,max_grids-1
     gg(n)%dxg = gg(n-1)%dxg*2.d0
     gg(n)%dyg = gg(n-1)%dyg*2.d0
#ifdef TWO_D
     gg(n)%kmax = 1
     gg(n)%dzg  = gg(n-1)%dzg
#else 
     gg(n)%dzg = gg(n-1)%dzg*2.d0
#endif
     ! take care of reductions
     if( is_top( gg(n-1) ) ) then
        ! all done - clear rest of grids
        do ii=n,max_grids-1
           gg(ii)%xprocs = 0 
           gg(ii)%yprocs = 0 
           gg(ii)%zprocs = 0 
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
     else if (gg(n-1)%imax>min_psize .and. gg(n-1)%jmax>min_psize &
#ifndef TWO_D
          .and. gg(n-1)%kmax>min_psize &
#endif
        ) then
        !     normal reduction, no split yet	      
        gg(n)%imax=gg(n-1)%imax/2
        gg(n)%jmax=gg(n-1)%jmax/2
#ifndef TWO_D
        gg(n)%kmax=gg(n-1)%kmax/2
#endif
        gg(n)%comm3d = comm3d
        gg(n)%comm = MPI_COMM_WORLD
        gg(n)%iprocx = gg(n-1)%iprocx
        gg(n)%iprocy = gg(n-1)%iprocy
        gg(n)%iprocz = gg(n-1)%iprocz
        gg(n)%xprocs = gg(n-1)%xprocs
        gg(n)%yprocs = gg(n-1)%yprocs
        gg(n)%zprocs = gg(n-1)%zprocs
         
        gg(n)%left = gg(n-1)%left
        gg(n)%right = gg(n-1)%right
        gg(n)%top = gg(n-1)%top
        gg(n)%bottom = gg(n-1)%bottom
        gg(n)%forward = gg(n-1)%forward
        gg(n)%behind = gg(n-1)%behind
     
        gg(n)%loc_comm = mpi_comm_null
     else  ! gg(n-1)%imax == min_psize ...
        !     reduce number of procs on grid by 2
        gg(n)%xprocs = gg(n-1)%xprocs/2
        gg(n)%yprocs = gg(n-1)%yprocs/2
        gg(n)%iprocx = (gg(n-1)%iprocx-1)/2 + 1
        gg(n)%iprocy = (gg(n-1)%iprocy-1)/2 + 1
#ifdef TWO_D
        gg(n)%zprocs = gg(n-1)%zprocs
        gg(n)%iprocz = gg(n-1)%iprocz
#else 
        gg(n)%zprocs = gg(n-1)%zprocs/2
        gg(n)%iprocz = (gg(n-1)%iprocz-1)/2 + 1
#endif
        if(gg(n)%zprocs==0 .or. gg(n)%yprocs==0 .or. &
             gg(n)%xprocs==0 ) then
           gg(n)%zprocs=0; gg(n)%yprocs=0; gg(n)%xprocs=0
           if(mype==0) write(6,*)'[',mype,'] domain is too thin'
           exit ! domain is too thin
        endif
        ! zero based proc ID
        proc_id = (gg(n)%iprocx-1)*gg(n)%yprocs*gg(n)%zprocs &
             + (gg(n)%iprocy-1)*gg(n)%zprocs + (gg(n)%iprocz-1)
        ii = mod(gg(n-1)%iprocx-1,2)
        jj = mod(gg(n-1)%iprocy-1,2)
        kk = mod(gg(n-1)%iprocz-1,2)
#ifdef TWO_D
        ii = jj + 2*ii    ! local zero based ID
#else 
        ii = kk + 2*jj + 4*ii ! local zero based ID
#endif
        call MPI_COMM_SPLIT(gg(n-1)%comm,ii,proc_id,gg(n)%comm,ierr)
        call MPI_COMM_SPLIT(gg(n-1)%comm,proc_id,ii,gg(n)%loc_comm,ierr)

        NProcAxis(1) = gg(n)%xprocs
        NProcAxis(2) = gg(n)%yprocs
        NProcAxis(3) = gg(n)%zprocs
        if(mype==0)write(0,*)'[',mype,']cart_create=',NProcAxis
        call MPI_cart_create( gg(n)%comm, ndims, NProcAxis, &
             periods, reorder, gg(n)%comm3D, ierr)
        if(ierr.ne.0) stop 'MPI_cart_create'
        
        ! debug              
        call MPI_comm_rank(gg(n)%comm3D,ii,ierr)
        call MPI_Cart_Coords(gg(n)%comm3D,ii,ndims,NProcAxis,ierr)
        if(mype==0)then
           write(0,*) '[',mype, '] cart rank',ii
           write(0,*) '[',mype, '] X',gg(n)%iprocx,NProcAxis(1)+1
           write(0,*) '[',mype, '] Y',gg(n)%iprocy,NProcAxis(2)+1
           write(0,*) '[',mype, '] Z',gg(n)%iprocz,NProcAxis(3)+1
           ii = mod(gg(n-1)%iprocx-1,2)
           jj = mod(gg(n-1)%iprocy-1,2)
           write(0,*) '[',mype,'] ii=',ii,'jj=',jj,'proc_id=',proc_id
           write(0,*) '[',mype, '] nx', gg(n)%xprocs
           write(0,*) '[',mype, '] ny', gg(n)%yprocs
           write(0,*) '[',mype, '] nz', gg(n)%zprocs
        endif        
        if(gg(n)%iprocx.ne.NProcAxis(1)+1) stop '%iprocx'
        if(gg(n)%iprocy.ne.NProcAxis(2)+1) stop '%iprocy'
        if(gg(n)%iprocz.ne.NProcAxis(3)+1) stop '%iprocz'
        
        !      write(0,*)'[',mype,']ipx=',gg(n)%iprocx,gg(n)%iprocy,gg(n)%iprocz
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
!!$  ixlo=1-nghost
!!$  iylo=1-nghost
!!$  izlo=1-nghost
!!$  ixhi=nxlocal+nghost
!!$  iyhi=nylocal+nghost
!!$  izhi=nzlocal+nghost 
!!$#ifdef TWO_D
!!$  izlo=1; izhi=nzlocal ! ???
!!$#endif
  !     Allocate mesh
!!$  allocate(xc(ixlo:ixhi))
!!$  allocate(yc(iylo:iyhi))
!!$  allocate(zc(izlo:izhi))

  enddo
end subroutine new_grids_private
