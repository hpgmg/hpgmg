!-----------------------------------------------------------------
!     Ravi Samtaney & Mark Adams
!     Copyright 2014
!-----------------------------------------------------------------
subroutine solve(g,Apply1,Apply2,Relax1,Res0,errors,nViters)
  use GridModule
  use pms, only:dom_max_grids,nvcycles,nfcycles,rtol,verbose,dom_max_sr_grids,nsr,ncgrids
  use mpistuff
  use error_data_module
  implicit none
  !
  type(pe_patch),intent(in)::g(-dom_max_sr_grids:dom_max_grids-1)
  type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
  external Apply1, Apply2, Relax1
  double precision,intent(out)::Res0
  integer,intent(out):: nViters
  !     Local Vars
  integer:: iter
  double precision:: Res,Reslast
  double precision,dimension(&
       g(0)%ilo:g(0)%ihi, g(0)%jlo:g(0)%jhi,&
       g(0)%klo:g(0)%khi, nvar) :: L2u_f,rhs,ux
  double precision norm

  !call Apply2(L2u_f,ux,g(0)) ! ux==0 
  call formRHS(rhs,g(0))
  Res0 = norm(rhs,g(0),2)
  if(verbose.gt.1) then
     call formExactU(L2u_f,g(0))
     Reslast = norm(L2u_f,g(0),2) 
     if(mype==0) write(6,'(A,E14.6,A,E14.6,A,I5)') '0) solve: |f|_2=',&
          Res0,', |u_exact|_2=',Reslast,', nx=',g(0)%imax*g(0)%npex
  end if
  Reslast = Res0
  ! FMG (SR coarse grid) solve
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
     ux=0.d0
  end if
 
  ! finish with V-cycle - no V with SR
  nViters = 0
  do while(Res/Res0 > rtol .and. iter < nvcycles .and. nsr==0)
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
! MGF: recursive F-cycle with data on stack
!-----------------------------------------------------------------------
recursive subroutine MGF(ux,rhs,g,lev,Apply1,Apply2,Relax1,errors)
  use GridModule
  use mpistuff,only:mype
  use error_data_module
  use pms, only:dom_max_grids,nfmgvcycles,ncoarsesolveits,verbose,nvcycles,dom_max_sr_grids,nsr,ncgrids
  implicit none
  integer,intent(in) :: lev
  type(pe_patch),intent(in)::g(-dom_max_sr_grids:dom_max_grids-1)
  type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
  double precision,intent(in),dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: rhs
  double precision,intent(out),dimension(&
       g(lev)%ilo:g(lev)%ihi, g(lev)%jlo:g(lev)%jhi,& 
       g(lev)%klo:g(lev)%khi, nvar) :: ux
  external Apply1, Apply2, Relax1
  ! Local Vars - Coarse. This is new data on the stack
  integer :: ii,jj,kk
  double precision :: norm,tt
  ! Coarse data. This is new data on the stack
  double precision,dimension(g(lev+1)%ilo:g(lev+1)%ihi,g(lev+1)%jlo:g(lev+1)%jhi,&
       g(lev+1)%klo:g(lev+1)%khi, nvar):: uxC,FC
  double precision,dimension(g(lev)%ilo:g(lev)%ihi,g(lev)%jlo:g(lev)%jhi,&
       g(lev)%klo:g(lev)%khi, nvar)::tmp

  ux=0.d0 ! everthing uses initial solution except FMG
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
  end if

  ! form errors & convergance measure, recursive so goes from coarse to fine (good)
  call formExactU(tmp,g(lev))
  errors(err_lev)%uerror = norm(tmp-ux,g(lev),error_norm)
  call formGradError(tmp,ux,g(lev),errors(err_lev)%graduerr,error_norm)
  call Apply2(tmp,ux,g(lev))
  errors(err_lev)%resid = norm(tmp-rhs,g(lev),2)
 
  if (verbose.gt.0) then
     if(mype==0)write(6,'(A,I2,A,E14.6,A,E14.6,A,I8,I8,I8,A,I4,I4,I4)') &
          '     lev=',lev,') MGF |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror,&
          ', n=',g(lev)%imax*g(lev)%npex,g(lev)%jmax*g(lev)%npey,g(lev)%kmax*g(lev)%npez,&
          ', npe=',g(lev)%npex,g(lev)%npey,g(lev)%npez
  end if
  
  err_lev = err_lev + 1
  return
end subroutine MGF
!-----------------------------------------------------------------------
recursive subroutine MGV(ux,rhs,g,lev,Apply,Relax)
  use GridModule
  use mpistuff, only:mype
  use pms, only:dom_max_grids,verbose,ncoarsesolveits,ncycles,nsmoothsdown,&
       nsmoothsup,pe_min_sz,dom_max_sr_grids
  !
  implicit none
  integer,intent(in):: lev
  type(pe_patch),intent(in)::g(-dom_max_sr_grids:dom_max_grids-1)
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
     if (verbose.gt.2) then
        tt = norm(ux,g(lev),2); if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: coarse u_1=',tt
     end if
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
subroutine destroy_grids_private(gg)
  use GridModule
  use mpistuff
  use pms, only:pe_min_sz,dom_max_sr_grids,dom_max_grids,nsr
  implicit none
  type(pe_patch):: gg(-dom_max_sr_grids:dom_max_grids-1)
  ! 
  integer :: n

  do n=-nsr,dom_max_grids-1,1
     if (n.le.0) then
        ! nothing to delete for SR and grid 0
     else if (gg(n)%imax==0) then
        exit ! past top
     else if (gg(n-1)%imax > pe_min_sz .and. gg(n-1)%jmax > pe_min_sz &
#ifdef TWO_D
          ) then
#else
        .and. gg(n-1)%kmax > pe_min_sz ) then
#endif
        ! normal reduction, no split yet -- nothing to destroy
     else 
        ! reduce number of pes on grid by 2
        if (gg(n)%npex .ne. gg(n-1)%npex ) then
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
subroutine new_grids_private(gg,NPeAxis,iPeAxis,nx,ny,nz,&
     nxlocal,nylocal,nzlocal,comm3d)
  use GridModule
  use mpistuff
  use pms
  use domain
  implicit none
  integer:: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
  integer:: NPeAxis(3),iPeAxis(3)
  type(pe_patch),intent(out):: gg(-dom_max_sr_grids:dom_max_grids-1)
  integer :: n,ndims,ii,jj,kk,pe_id,rank
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
  !       ipex ...	
  gg(0)%ipex = iPeAxis(1)
  gg(0)%ipey = iPeAxis(2)
  gg(0)%ipez = iPeAxis(3)
  
  gg(0)%npex = NPeAxis(1) 
  gg(0)%npey = NPeAxis(2) 
  gg(0)%npez = NPeAxis(3) 
  ! 
  gg(0)%imax=nxlocal
  gg(0)%jmax=nylocal
#ifdef TWO_D
  gg(0)%kmax=1
#else 
  gg(0)%kmax=nzlocal
#endif
  ! start at 1 as 0 was just done
  ncgrids = 1
  do n=1,dom_max_grids-1
     gg(n)%dxg = gg(n-1)%dxg*mg_ref_ratio
     gg(n)%dyg = gg(n-1)%dyg*mg_ref_ratio
#ifdef TWO_D
     gg(n)%kmax = 1
     gg(n)%dzg  = gg(n-1)%dzg
#else 
     gg(n)%dzg = gg(n-1)%dzg*mg_ref_ratio
#endif
     ! take care of reductions
     if( is_top(gg(n-1)) ) then
        ! all done - clear rest of grids
        do ii=n,dom_max_grids-1
           gg(ii)%npex = 0 
           gg(ii)%npey = 0 
           gg(ii)%npez = 0 
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
     else if (gg(n-1)%imax>pe_min_sz .and. gg(n-1)%jmax>pe_min_sz &
#ifndef TWO_D
          .and. gg(n-1)%kmax>pe_min_sz &
#endif
          ) then
        !     normal reduction, no split yet	      
        gg(n)%imax=gg(n-1)%imax/mg_ref_ratio
        gg(n)%jmax=gg(n-1)%jmax/mg_ref_ratio
#ifndef TWO_D
        gg(n)%kmax=gg(n-1)%kmax/mg_ref_ratio
#else
        gg(n)%kmax=gg(n-1)%kmax ! 1
#endif
        gg(n)%comm3d = comm3d
        gg(n)%comm = MPI_COMM_WORLD
        gg(n)%ipex = gg(n-1)%ipex
        gg(n)%ipey = gg(n-1)%ipey
        gg(n)%ipez = gg(n-1)%ipez
        gg(n)%npex = gg(n-1)%npex
        gg(n)%npey = gg(n-1)%npey
        gg(n)%npez = gg(n-1)%npez
        
        gg(n)%left = gg(n-1)%left
        gg(n)%right = gg(n-1)%right
        gg(n)%top = gg(n-1)%top
        gg(n)%bottom = gg(n-1)%bottom
        gg(n)%forward = gg(n-1)%forward
        gg(n)%behind = gg(n-1)%behind
        gg(n)%loc_comm = mpi_comm_null
     else  ! gg(n-1)%imax == min_psize ...
        !     reduce number of pes on grid by mg_ref_ratio
        if(gg(n-1)%npez==1 .or. gg(n-1)%npey==1 &
#ifndef TWO_D
             .or. gg(n-1)%npex==1 ) then ! one pe reduction, copy comm stuff
#else
           ) then
#endif
           gg(n)%npex = gg(n-1)%npex
           gg(n)%npey = gg(n-1)%npey
           gg(n)%ipex = gg(n-1)%ipex
           gg(n)%ipey = gg(n-1)%ipey
           gg(n)%npez = gg(n-1)%npez
           gg(n)%ipez = gg(n-1)%ipez
           
           gg(n)%imax=gg(n-1)%imax/mg_ref_ratio
           gg(n)%jmax=gg(n-1)%jmax/mg_ref_ratio
#ifndef TWO_D
           gg(n)%kmax=gg(n-1)%kmax/mg_ref_ratio
#else
           gg(n)%kmax=1
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
              write(0,*) '[',mype, '] X:',gg(n)%ipex
              write(0,*) '[',mype, '] Y:',gg(n)%ipey
              write(0,*) '[',mype, '] Z:',gg(n)%ipez
              write(0,*) '[',mype, '] nx', gg(n)%npex
              write(0,*) '[',mype, '] ny', gg(n)%npey
              write(0,*) '[',mype, '] nz', gg(n)%npez
           endif
        else ! normal split
           gg(n)%npex = gg(n-1)%npex/mg_ref_ratio
           gg(n)%npey = gg(n-1)%npey/mg_ref_ratio
           gg(n)%ipex = (gg(n-1)%ipex-1)/mg_ref_ratio + 1
           gg(n)%ipey = (gg(n-1)%ipey-1)/mg_ref_ratio + 1
#ifdef TWO_D
           gg(n)%npez = gg(n-1)%npez
           gg(n)%ipez = gg(n-1)%ipez
#else
           gg(n)%npez = gg(n-1)%npez/mg_ref_ratio
           gg(n)%ipez = (gg(n-1)%ipez-1)/mg_ref_ratio + 1
#endif
           ! zero based pe ID
           pe_id = (gg(n)%ipex-1)*gg(n)%npey*gg(n)%npez &
                + (gg(n)%ipey-1)*gg(n)%npez + (gg(n)%ipez-1)
           ii = mod(gg(n-1)%ipex-1,2)
           jj = mod(gg(n-1)%ipey-1,2)
           kk = mod(gg(n-1)%ipez-1,2)
#ifdef TWO_D
           ii = jj + 2*ii           ! local zero based ID
#else 
           ii = kk + 2*jj + 4*ii    ! local zero based ID
#endif
           call MPI_COMM_SPLIT(gg(n-1)%comm,ii,pe_id,gg(n)%comm,ierr)
           call MPI_COMM_SPLIT(gg(n-1)%comm,pe_id,ii,gg(n)%loc_comm,ierr)

           NPeAxis(1) = gg(n)%npex
           NPeAxis(2) = gg(n)%npey
           NPeAxis(3) = gg(n)%npez
           
           call MPI_cart_create( gg(n)%comm, ndims, NPeAxis, &
                periodic, .false., gg(n)%comm3D, ierr)
           
           ! debug              
           call MPI_comm_rank(gg(n)%comm3D,ii,ierr)
           call MPI_Cart_Coords(gg(n)%comm3D,ii,ndims,NPeAxis,ierr)
           if(gg(n)%ipex.ne.NPeAxis(1)+1) stop '%ipex'
           if(gg(n)%ipey.ne.NPeAxis(2)+1) stop '%ipey'
           if(gg(n)%ipez.ne.NPeAxis(3)+1) stop '%ipez'
           if (mype==npe/2.and.verbose.gt.4)then
              write(0,*) '[',mype, '] cart rank',ii
              write(0,*) '[',mype, '] X:',gg(n)%ipex,NPeAxis(1)+1
              write(0,*) '[',mype, '] Y:',gg(n)%ipey,NPeAxis(2)+1
              write(0,*) '[',mype, '] Z:',gg(n)%ipez,NPeAxis(3)+1
              ii = mod(gg(n-1)%ipex-1,2)
              jj = mod(gg(n-1)%ipey-1,2)
              write(0,*) '[',mype,'] ii=',ii,'jj=',jj,'pe_id=',pe_id
              write(0,*) '[',mype, '] nx', gg(n)%npex
              write(0,*) '[',mype, '] ny', gg(n)%npey
              write(0,*) '[',mype, '] nz', gg(n)%npez
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
     if(gg(n)%npez==0 .or. gg(n)%npey==0 .or. &
          gg(n)%npex==0 ) then
        gg(n)%npez=0; gg(n)%npey=0; gg(n)%npex=0
        if(mype==0) write(6,*)'[',mype,'] domain is too thin',gg(n)%npex==0,gg(n)%npey==0,gg(n)%npez==0
        stop 'domain is too thin' ! domain is too thin
     endif
     
     ! print topology
     if(mype==0.and.verbose.gt.1) then
        write(6,*) '[',mype,'] level ',n,', nxl=',gg(n)%imax,'npx=',gg(n)%npex
        if (mype==-1) then
           write(6,*) '[',mype,'] ipx=',gg(n)%ipex, ',ipy=',gg(n)%ipey,',ipz= ',gg(n)%ipez
           write(6,*) '[',mype,'] npx=',gg(n)%npex, ',npy=',gg(n)%npey,',npz=',gg(n)%npez
           write(6,*) '[',mype,'] nxl=',gg(n)%imax,   ',nyl=',gg(n)%jmax,  ',nzl=',gg(n)%kmax
        endif
     end if
     ncgrids = ncgrids + 1 ! keep track for ease
  enddo ! non- finest grid construction

  ! SR stuff, from coarse to fine
  nsr = 0 ! global var set here
  do n=-1,-dom_max_sr_grids,-1     
     ! set everthing
     gg(n)%loc_comm = mpi_comm_null
     gg(n)%comm3d = mpi_comm_null
     gg(n)%comm = mpi_comm_null
     gg(n)%ipex = 0
     gg(n)%ipey = 0
     gg(n)%ipez = 0
     gg(n)%npex = 0
     gg(n)%npey = 0
     gg(n)%npez = 0
     
     gg(n)%imax=0 ! used as flag

     ! size of grid
     nxlocal = nxlocal*mg_ref_ratio 
     nylocal = nylocal*mg_ref_ratio
     nzlocal = nzlocal*mg_ref_ratio

     ! are we going to make an SR grid?
     if ( nxlocal .le. sr_max_loc_sz .or. nylocal .le. sr_max_loc_sz &
#ifdef TWO_D
          ) then
#else 
        .or. nzlocal .le. sr_max_loc_sz) then
#endif
        nsr = nsr + 1 ! have sr
        ! dx
        gg(n)%dxg = gg(n-1)%dxg*mg_ref_ratio
        gg(n)%dyg = gg(n-1)%dyg*mg_ref_ratio
#ifdef TWO_D
        gg(n)%kmax = 1
        gg(n)%dzg  = gg(n-1)%dzg
#else 
        gg(n)%dzg = gg(n-1)%dzg*mg_ref_ratio
#endif
        ! new (valid) size
        gg(n)%imax=nxlocal
        gg(n)%jmax=nylocal
#ifdef TWO_D
        gg(n)%kmax=1
#else 
        gg(n)%kmax=nzlocal
#endif
     else
        exit
     end if
  end do

  ! add buffer/ghosts, the size
  ii = sr_base_bufsz - sr_bufsz_inc ! the buffer schedual
  do n=-nsr,ncgrids-1,1
     if (n.lt.0) then
        ii = ii + sr_bufsz_inc ! number of buffer cells
        gg(n)%ivallo= ii+1
        gg(n)%jvallo=ii+1
        gg(n)%ivalhi=gg(n)%imax+ii
        gg(n)%jvalhi=gg(n)%jmax+ii
        gg(n)%imax=gg(n)%imax+2*ii ! new comp area
        gg(n)%jmax=gg(n)%jmax+2*ii
#ifndef TWO_D
        gg(n)%kvallo=ii+1
        gg(n)%kvalhi=gg(n)%kmax+ii
        gg(n)%kmax=gg(n)%kmax+2*ii
#else
        gg(n)%kmax=1
        gg(n)%kvallo=1
        gg(n)%kvalhi=1
#endif 
     else ! normal grid valid region
        gg(n)%ivallo=1
        Gg(n)%jvallo=1
        gg(n)%kvallo=1
        gg(n)%ivalhi=gg(n)%imax
        Gg(n)%jvalhi=gg(n)%jmax
        gg(n)%kvalhi=gg(n)%kmax
        ii = 0 ! no more SR buffs (for print)
     end if
     ! data size
     gg(n)%ilo=-nsg+1
     gg(n)%ihi= gg(n)%imax+nsg
     gg(n)%jlo=-nsg+1
     gg(n)%jhi=gg(n)%jmax+nsg
#ifdef TWO_D
     gg(n)%klo=1
     gg(n)%khi=1
#else 
     gg(n)%klo=-nsg+1
     gg(n)%khi=gg(n)%kmax+nsg
#endif
     if (verbose.gt.0) then
        if (mype==0) write (6,'(A,I2,A,I2,I2,I2,A,I4,I4,I4,A,I4,I4,I4,A,I3)'),&
             'new_grids:',n,') valid lo:',gg(n)%ivallo,gg(n)%jvallo,gg(n)%kvallo&
             ,' valid hi:',gg(n)%ivalhi,gg(n)%jvalhi,gg(n)%kvalhi&
             ,' max:',gg(n)%imax,gg(n)%jmax,gg(n)%kmax,' sr nbuf',ii
     end if
  end do

end subroutine new_grids_private
   
