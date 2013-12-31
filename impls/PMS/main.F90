!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
program PMS
  use GridModule, only:reorder
  use mpistuff
  implicit none
  integer:: rank,comm3d
  integer:: ndims,NProcAxis(3),iProcAxis(3),xprocs,yprocs,zprocs,nx,ny,nz,nxlocal,nylocal,nzlocal
  logical:: periods(3)
  ! MPI init - globals
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, mype, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, npe, ierr)

  call InitIO()

  ! profiling package init
#ifdef V_T
  include 'VT.inc'
  call VTTRACEOFF(ierr)
#elif PMS_HAVE_PAPI
  include 'f90papi.h'
  INTEGER retval
  call PAPIf_num_counters( ncounters )
  if ( ncounters.le.0) then
     print *, '   error: PAPIf_library_init',PAPI_VER_CURRENT
     call abort
  end if
  !
  call PAPIf_query_event(PAPI_FP_OPS, retval)
  if (retval .NE. PAPI_OK) then
     print *, '       error: PAPIf_query_event: ',retval,PAPI_OK
     call abort
  end if
  ! create the eventset
  fp_event = PAPI_NULL
  call PAPIF_create_eventset(fp_event, retval)
  if (retval .ne. PAPI_OK) then
     print *, 'Error in subroutine PAPIF_create_eventset'
     call abort
  end if
  ! set counters to count flops
  call PAPIF_add_event(fp_event, PAPI_FP_OPS, retval)
  if (retval .NE. PAPI_OK) then
     print *, "Abort After PAPIF_add_events: ", retval
     call abort
  endif
#endif

  ! get fine grid parameters and problem type
  call get_params(npe, xprocs, yprocs, zprocs, nx, ny, nz, &
       nxlocal, nylocal, nzlocal)

  ! Initialize the fine grid communicator
  ndims = 3
  NProcAxis(1) = xprocs
  NProcAxis(2) = yprocs
  NProcAxis(3) = zprocs
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
  call MPI_cart_create(MPI_COMM_WORLD,ndims,NProcAxis,periods,&
       reorder,comm3d,ierr)
  call MPI_comm_rank(comm3d,rank,ierr) ! new rank if reorder==true -- OK!!!

  if(mype==0) then
     write(6,*) 'main: mype=',mype,', comm3D rank = ',rank
  endif

  call MPI_Cart_Coords(comm3d,rank,ndims,iProcAxis,ierr)
  iProcAxis(1) = iProcAxis(1) + 1 ! one based
  iProcAxis(2) = iProcAxis(2) + 1
  iProcAxis(3) = iProcAxis(3) + 1

  ! do it
  call go(NProcAxis, iProcAxis, nx, ny, nz, nxlocal, nylocal, nzlocal, comm3d)

  ! end it
  call FinalizeIO()
  call MPI_Finalize(ierr)
end program PMS
!-----------------------------------------------------------------
!  go - sets up grids and calls the main 'driver'
!-----------------------------------------------------------------
subroutine go(NProcAxis,iProcAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)
  use GridModule
  implicit none
  integer,intent(in):: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
  integer :: NProcAxis(3),iProcAxis(3)
  !
  type(proc_patch)::grids(0:max_grids-1)
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

  call SetupDomain() ! sets size of domain
  
  call new_grids(grids,NProcAxis,iProcAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)

  ! do it again
  call driver(grids)

  call destroy_grids_private(grids,max_grids)
  
end subroutine go
!-----------------------------------------------------------------
!  driver  
!-----------------------------------------------------------------
subroutine driver(grids)
  !
  use iounits
  use mpistuff, only:mype
  use domain, only:problemType
  use GridModule
  use error_data_module
  implicit none
  !
  type(proc_patch)::grids(0:max_grids-1)
  !
  integer:: isolve,maxiter,i,coarsest_grid,nVcycles
  !integer:: output_flag, binary_flag
  integer,parameter:: dummy_size=1
  double precision:: dummy(dummy_size),res0,rateu,rategrad
  type(error_data)::errors(0:max_grids-1+nfvcycles)
  external Apply_const_Lap
  external GSRB_const_Lap
  double precision,parameter:: log2r = 1.d0/log(2.d0);
print *,1111
  !output_flag=0 ! parameters
  !binary_flag=1
  maxiter = 1   ! parameter (10-1000)

  do isolve=1,maxiter
     ! flush cache
     do i=1,dummy_size
        dummy(i) = 1.235d3
     end do

     select case (problemType)
     case(0)
  
        call solve(grids,Apply_const_Lap,Apply_const_Lap,GSRB_const_Lap,res0,errors,nVcycles)

     case(1)

        if(mype==0) write(6,*) 'problemType 1 not implemented -- todo'

     end select

     if(mype==0) then
        write(6,*) isolve,') error=',errors(0)%uerror,', rel resid=',&
             errors(0)%resid/res0,', grad(u) error=',errors(0)%graduerr
        write(ihis,*) isolve,errors(0)%uerror,errors(0)%graduerr,errors(0)%resid/res0
        coarsest_grid = -1
        do i=max_grids-1,1,-1
           ! go from coarse to fine.  Need to pass through top empty 'coarse' grids
           if( .not. is_top(grids(i)) ) then
              rateu = errors(i)%uerror/errors(i+1)%uerror
              rategrad = errors(i)%graduerr/errors(i+1)%graduerr
              write(6,'(I5,A8,I2,A20,E14.6,A10,E14.6)') isolve,': level ',i, ': converg order u=', &
                   log(rateu)*log2r,', grad(u)=',log(rategrad)*log2r
              if (i == 1) then ! print fine grid 
                 write(iconv,*) isolve,log(rateu)*log2r,log(rategrad)*log2r,errors(i)%uerror,&
                      errors(i)%graduerr,errors(i)%resid
              end if
              if (coarsest_grid == -1) coarsest_grid = i+1
              print *,'driver: do level i=',i
           end if
        end do
        print *, 'driver: coarsest grid = ',coarsest_grid
        do i=coarsest_grid,coarsest_grid+nVcycles
           print *,'driver: do V-cycles i=',i
           write(6,'(A13,I2,A20,E14.6,A8,E14.6)') 'driver: level', i, ': error=',errors(i)%uerror,&
                ': resid=',errors(i)%resid
           write(iconv,*) isolve,0.d0,0.d0,errors(i)%uerror,errors(i)%graduerr,errors(i)%resid
        end do
     endif
  enddo

!!$  if(output_flag.eq.2) then
!!$     if(binary_flag.nq.0) then
!!$        if(mype==0) then
!!$           write(6,*) 'Binary output at',timeStep, ttot
!!$        endif
!!$        Call WriteSiloBinaryFileParallel(ux,timeStep)
!!$     endif
!!$  endif
  
  return
end subroutine driver
!-----------------------------------------------------------------
! get command line args, set problemType and fine grid sizes
!-----------------------------------------------------------------       
subroutine get_params(nprocs, xprocs, yprocs, zprocs, &
     nx, ny, nz, nxlocal, nylocal, nzlocal)
  use domain, only:problemType
  use mpistuff
#ifndef HAVE_COMM_ARGS
  use f2kcli        ! command line args class
#endif
  !use mpistuff, only:mype
  implicit none
  integer,intent(out) :: xprocs,yprocs,zprocs,nx,ny,nz,nxlocal,nylocal,nzlocal
  integer,intent(in) :: nprocs
  !       
  CHARACTER(LEN=256) :: LINE
  CHARACTER(LEN=10)  :: CMD,CMD2
  INTEGER            :: NARG,IARG,fact
  integer*8 :: ti
  double precision t2,t1,iarray(12)
  !       default input args: ?local, n?, ?procs
  xprocs=-1
  yprocs=-1
  zprocs=-1
  nx=-1
  ny=-1
  nz=-1
  nxlocal=-1
  nylocal=-1
  nzlocal=-1
  !       Problem type 
  !       1: 
  problemType = 1
  if( mype==0 ) then
     NARG = COMMAND_ARGUMENT_COUNT() !
     DO IARG = 1,NARG,2
        CALL GET_COMMAND_ARGUMENT(IARG,CMD)
        CALL GET_COMMAND_ARGUMENT(IARG+1,CMD2)
        if( CMD == '-nx' ) READ (CMD2, '(I10)') nx
        if( CMD == '-ny' ) READ (CMD2, '(I10)') ny
        if( CMD == '-nz' ) READ (CMD2, '(I10)') nz
        
        if( CMD == '-nxlocal' ) READ (CMD2, '(I10)') nxlocal
        if( CMD == '-nylocal' ) READ (CMD2, '(I10)') nylocal
        if( CMD == '-nzlocal' ) READ (CMD2, '(I10)') nzlocal
        
        if( CMD == '-xprocs' ) READ (CMD2, '(I10)') xprocs
        if( CMD == '-yprocs' ) READ (CMD2, '(I10)') yprocs
        if( CMD == '-zprocs' ) READ (CMD2, '(I10)') zprocs
        
        if( CMD == '-problem' ) READ (CMD2, '(I10)') problemType
     END DO
#ifdef TWO_D
     nz = 1; nzlocal = 1; zprocs = 1;
#endif
     iarray(1) = nx
     iarray(2) = ny
     iarray(3) = nz
     
     iarray(4) = nxlocal
     iarray(5) = nylocal
     iarray(6) = nzlocal
     
     iarray(7) = xprocs
     iarray(8) = yprocs
     iarray(9) = zprocs
     
     iarray(10) = problemType
  
     call MPI_Bcast(iarray,10,MPI_DOUBLE_PRECISION, 0, &
          mpi_comm_world, ierr)
  else
     call MPI_Bcast(iarray,10,MPI_DOUBLE_PRECISION, 0, &
          mpi_comm_world, ierr)
     nx = int(iarray(1))
     ny = int(iarray(2))
     nz = int(iarray(3))
     
     nxlocal = int(iarray(4))
     nylocal = int(iarray(5))
     nzlocal = int(iarray(6))
     
     xprocs = int(iarray(7))
     yprocs = int(iarray(8))
     zprocs = int(iarray(9))
     
     problemType = int(iarray(10))
  endif
  ! get ijk procs
  t1=nprocs
#ifdef TWO_D
  zprocs = 1
  nz = 1
  nzlocal = ????
  if (xprocs == -1) then
     if( abs(dsqrt(t1) - floor(dsqrt(t1))) > 0.d0 ) then
        if(abs(dsqrt(t1/2.d0)-floor(dsqrt(t1/2.d0)))>0.d0)then
           if(abs(dsqrt(t1/4.d0)-floor(dsqrt(t1/4.d0)))>0.d0)then
              fact=8    ! hope this works
           else
              fact=4
           endif
        else
           fact=2
        endif
     else
        fact=1
     endif
     xprocs=fact*floor(dsqrt(t1/fact))
     yprocs=t1/xprocs
     if (nprocs .ne. zprocs*xprocs*yprocs) then
        write(6,*) 'INCORRECT NP:', NPROCS,'xprocs=',xprocs,&
             'yprocs=',yprocs,'zprocs=',zprocs,'nprocs=',nprocs
        stop
     endif
  else
     yprocs=t1/xprocs
     if(yprocs<1)stop'too many xprocs???'
     fact = xprocs/yprocs
     if( fact < 1 ) then
        write(6,*) 'xprocs must be greater than yprocs', NPROCS,&
       'xprocs=',xprocs,&
       'yprocs=',yprocs,'zprocs=',zprocs,'nprocs=',nprocs
        stop
     endif
  endif
#else
  if (xprocs == -1) then
     t1 = nprocs
     if(nprocs==1) then
        xprocs = 1; yprocs = 1; zprocs = 1;
     elseif(nx==ny .and. ny==nz)then
        xprocs = floor(t1**.3333333334d0)
        yprocs = xprocs
        zprocs = nprocs/(xprocs*yprocs)
     elseif(nz==ny)then
        t1 = ny*nprocs; t2 = nx; t1 = t1/t2 
        yprocs = floor(t1**.333333334d0)
        zprocs = yprocs
        xprocs = nprocs/(yprocs*zprocs) 
     elseif(nz==nx)then
        t1 = nz*nprocs; t2 = ny; t1 = t1/t2 
        xprocs = floor(t1**.333333334d0)
        zprocs = xprocs
        yprocs = nprocs/(xprocs*zprocs) 
     elseif(nx==ny)then
        t1 = ny*nprocs; t2 = nz; t1 = t1/t2 
        yprocs = floor(t1**.333333334d0)
        xprocs = yprocs
        zprocs = nprocs/(yprocs*xprocs) 
     else
        stop 'need to specify x/y/z proc'
     endif
  endif
#endif     
  if (nprocs .ne. zprocs*xprocs*yprocs) then
     write(6,*) 'INCORRECT NP: xprocs=',xprocs, &
          'yprocs=',yprocs,'zprocs=',zprocs,'nprocs=',nprocs
     stop
  endif
  ! x
  if (nxlocal.eq.-1 .and. nx.eq.-1) nxlocal = 32        ! default
  if (nxlocal.eq.-1 .and. nx.ne.-1) nxlocal = nx/xprocs 
  if (nxlocal.ne.-1 .and. nx.eq.-1) nx = nxlocal*xprocs 
  if (nx .ne. nxlocal*xprocs) stop 'nx .ne. nxlocal*xprocs'
  ! y
  if (nylocal.eq.-1 .and. ny.eq.-1) nylocal = nxlocal   ! default square
  if (nylocal.eq.-1 .and. ny.ne.-1) nylocal = ny/yprocs 
  if (nylocal.ne.-1 .and. ny.eq.-1) ny = nylocal*yprocs 
  if (ny .ne. nylocal*yprocs) stop 'ny .ne. nylocal*yprocs'
  ! z
  if (nzlocal.eq.-1 .and. nz.eq.-1) nzlocal = nxlocal   ! default square
  if (nzlocal.eq.-1 .and. nz.ne.-1) nzlocal = nz/zprocs 
  if (nzlocal.ne.-1 .and. nz.eq.-1) nz = nzlocal*zprocs 
  if (nz .ne. nzlocal*zprocs) stop 'nz .ne. nzlocal*zprocs'
  ! print topology
  if(mype==0) then
     write(6,*) '[',mype,'] nx=',nx,',ny=',ny,',nz=',nz
     write(6,*)'[',mype,'] npx=',xprocs,',npy=',yprocs,',npz=',zprocs
     write(6,*)'[',mype,'] nxl=',nxlocal,',nyl=',nylocal,',nzl=',nzlocal
  endif
end subroutine get_params
!-----------------------------------------------------------------------
subroutine formExactU(exact,g)
  use GridModule
  use domain
  implicit none
  type(proc_patch),intent(in):: g
  double precision,intent(out)::exact(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar) 
  !
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)

  ig = iglobal(g) ! one bases
  jg = jglobal(g)
  kg = kglobal(g)
  select case (problemType)
  case(0)
     do kk=1,g%kmax
        do jj=1,g%jmax
           do ii=1,g%imax
              coord(1) = xl+(ig+ii-1)*g%dxg-0.5*g%dxg
              coord(2) = yl+(jg+jj-1)*g%dyg-0.5*g%dyg
#ifndef TWO_D
              coord(3) = zl+(kg+kk-1)*g%dzg-0.5*g%dzg
#endif
              x2 = coord*coord
              a = x2*(x2-1.);              
#ifdef TWO_D
              exact(ii,jj,kk,1) = a(1)*a(2)
#else
              exact(ii,jj,kk,1) = a(1)*a(2)*a(3)
#endif
           end do
        end do
     end do
  case(1)
     exact = 0.d0  
     print *,'**************** formExact not done for proble type 1: todo'
  end select
end subroutine formExactU
!-----------------------------------------------------------------------
subroutine formGradError(ux,exact,g,gerr,order)
  use GridModule
  use domain
  implicit none
  type(proc_patch),intent(in):: g
  double precision,intent(in)::exact(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar) 
  double precision,intent(in)::ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar) 
  double precision,intent(out)::gerr
  integer,intent(in):: order
  !
  integer:: ii,jj,kk
  double precision::grad1(3),grad2(3)
  double precision::gradError(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  double precision norm

  select case (problemType)
  case(0)
     do kk=1,g%kmax
        do jj=1,g%jmax
           do ii=1,g%imax
              grad1(1) = 0.5d0*(ux(ii+1,jj,kk,1)-ux(ii-1,jj,kk,1))/g%dxg
              grad1(2) = 0.5d0*(ux(ii,jj+1,kk,1)-ux(ii,jj-1,kk,1))/g%dyg
              grad2(1) = 0.5d0*(exact(ii+1,jj,kk,1)-exact(ii-1,jj,kk,1))/g%dxg
              grad2(2) = 0.5d0*(exact(ii,jj+1,kk,1)-exact(ii,jj-1,kk,1))/g%dyg
              
#ifndef TWO_D
              grad1(3) = 0.5d0*(ux(ii,jj,kk+1,1)-ux(ii,jj,kk-1,1))/g%dzg
              grad2(3) = 0.5d0*(exact(ii,jj,kk+1,1)-exact(ii,jj,kk-1,1))/g%dzg
              gradError(ii,jj,kk,1) = sqrt((grad1(1)-grad2(1))**2 + (grad1(2)-grad2(2))**2 &
                   + (grad1(3)-grad2(3))**2)
#else
              gradError(ii,jj,kk,1) = sqrt((grad1(1)-grad2(1))**2 + (grad1(2)-grad2(2))**2)
#endif
           end do
        end do
     end do
     gerr = norm(gradError,g,order)
  case(1)
     gerr = 0.d0  
     print *,'**************** formGradError not done for proble type 1: todo'
  end select
end subroutine formGradError
!-----------------------------------------------------------------------
subroutine FormRHS(rhs,g)
  use GridModule
  use domain
  implicit none
  type(proc_patch),intent(in):: g
  double precision,intent(out)::rhs(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)
 
  ig = iglobal(g)
  jg = jglobal(g)
  kg = kglobal(g)
!!$  !     
!!$  do k=IZLO, IZHI,1
!!$     zc(k)=zl+(kg+k-1)*dz-0.5*dz
!!$  enddo
!!$  do j=IYLO, IYHI,1
!!$     yc(j)=yl+(jg+j-1)*dy-0.5*dy
!!$  enddo
!!$  do i=IXLO, IXHI,1
!!$     xc(i)=xl+(ig+i-1)*dx-0.5*dx
!!$  enddo
  rhs = 1.23456789d300 ! debug
  select case (problemType)
  case(0)
     ! x^4 - x^2
#ifdef XPERIODIC
     stop 'x^4 - x^2 XPERIODIC not defined'
#endif
#ifdef YPERIODIC
     stop 'x^4 - x^2 YPERIODIC not defined'
#endif
#ifdef ZPERIODICtod
     stop 'x^4 - x^2 ZPERIODIC not defined'
#endif
     ! RHS = Lap(x^4-x^2)
     do kk=1,g%kmax
        do jj=1,g%jmax
           do ii=1,g%imax
              coord(1) = xl+(ig+ii-1)*g%dxg-0.5*g%dxg
              coord(2) = yl+(jg+jj-1)*g%dyg-0.5*g%dyg
#ifdef TWO_D
              coord(3) = zl+(kg+kk-1)*g%dzg-0.5*g%dzg
#endif
              print *, "coord=",coord
              x2 = coord*coord
              print *, "x2=",x2
              a = x2*(x2-1.);              
              b = 12.*x2-2.;
              print *, "b=",b
#ifdef TWO_D
              rhs(ii,jj,kk,1) = b(1)*a(2) + a(1)*b(2)
#else
              rhs(ii,jj,kk,1) = b(1)*a(2)*a(3) + a(1)*b(2)*a(3) + a(1)*a(2)*b(3)
#endif
           end do
        end do
     end do
  case(1)
     !     for bubble
#ifndef XPERIODIC
     stop 'bubble not XPERIODIC not defined'
#endif
#ifndef YPERIODIC
     stop 'bubble not YPERIODIC defined'
#endif
#ifndef ZPERIODIC
     stop 'bubble not ZPERIODIC defined'
#endif
     ! RHS = bubble
     

     stop 'FormRHS not impl -- todo'

  end select
  return
end subroutine FormRHS
!-----------------------------------------------------------------
subroutine SetupDomain()
  use domain
  implicit none
  !
  !     Problem type 
  !     0: (x^4-x^2)
  !     1: bubble
  !
  if(problemType.eq.0) then
     xl=0.D0
     xr=1.D0
     yl=0.D0
     yr=1.D0
     zl=0.D0
     zr=1.D0
  else if(problemType.eq.1) then
     xl=0.D0
     xr=1.D0
     yl=0.D0
     yr=1.D0
     zl=0.D0
     zr=1.D0
  else
     ! first place it is used
     stop 'SetupDomain: invalid problem type'
  endif
  return
end subroutine SetupDomain
!-----------------------------------------------------------------------
!     InitIO
!-----------------------------------------------------------------------
subroutine InitIO()
  use iounits
  use mpistuff, only:mype
  implicit none
  if(mype==0) then
     open(ihis,file='Run.history', form='formatted')
     open(iconv,file='Convergence.history', form='formatted')
  endif
  return
end subroutine InitIO
!     
!-----------------------------------------------------------------------
subroutine FinalizeIO()
  use iounits
  use mpistuff, only:mype
  implicit none
  ! Close all open files
  if(mype==0) then
     close(ihis)
     close(iconv)
  endif
  return
end subroutine FinalizeIO
