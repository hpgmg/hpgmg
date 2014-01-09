!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
program PMS
  use domain, only:verbose
  use mpistuff
  implicit none
  integer:: rank,comm3d,ndims
  integer:: NProcAxis(3),iProcAxis(3),nprocx,nprocy,nprocz,nx,ny,nz,nxlocal,nylocal,nzlocal
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
  call get_params(npe, nprocx, nprocy, nprocz, nx, ny, nz, &
       nxlocal, nylocal, nzlocal)

  ! Initialize the fine grid communicator
  ndims = 3
  NProcAxis(1) = nprocx
  NProcAxis(2) = nprocy
  NProcAxis(3) = nprocz
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
       .true.,comm3d,ierr) ! reoder is ignored???
  call MPI_comm_rank(comm3d,rank,ierr) ! new rank?
  if((mype==0.and.verbose.gt.3).or.mype.ne.rank) then
     if (verbose.gt.0) write(6,*) 'main: warning old mype=',mype,', comm3D new rank = ',rank
  endif
  mype=rank

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
!  go - sets up grids and calls the main driver
!-----------------------------------------------------------------
subroutine go(NProcAxis,iProcAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)
  use GridModule
  use domain, only:dom_max_grids
  implicit none
  integer,intent(in):: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
  integer :: NProcAxis(3),iProcAxis(3)
  !
  type(proc_patch)::grids(0:dom_max_grids-1)

  call SetupDomain() ! sets size of domain

  call new_grids(grids,NProcAxis,iProcAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)

  ! do it 
  call driver(grids)

  call destroy_grids_private(grids,dom_max_grids)
  
end subroutine go
!-----------------------------------------------------------------
!  driver  
!-----------------------------------------------------------------
subroutine driver(grids)
  use iounits
  use mpistuff 
  use domain
  use GridModule
  use error_data_module
  implicit none
  !
  type(proc_patch)::grids(0:dom_max_grids-1)
  !
  integer:: isolve,ii,coarsest_grid,nViters,kk
  !integer:: output_flag, binary_flag
  integer,parameter:: dummy_size=1024*1024 ! cache flush array size
  double precision:: dummy(dummy_size),res0,rateu,rategrad
  type(error_data)::errors(0:dom_max_grids-1+nvcycles)
  external Apply_const_Lap
  external GSRB_const_Lap
  double precision,parameter:: log2r = 1.d0/log(2.d0);
#ifdef HAVE_PETSC

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call PetscLogEventRegister('HPGMG solve', 0, events(1), ierr)
  call PetscLogEventRegister('HPGMG   smooth', 0, events(2), ierr)
  call PetscLogEventRegister('HPGMG   restrict', 0, events(3), ierr)
  call PetscLogEventRegister('HPGMG   apply', 0, events(4), ierr)
  call PetscLogEventRegister('HPGMG   prolong', 0, events(5), ierr)
  call PetscLogEventRegister('HPGMG   BC & exch', 0, events(6), ierr)
  call PetscLogEventRegister('HPGMG   form u,..', 0, events(7), ierr)
  call PetscLogEventRegister('HPGMG   other', 0, events(8), ierr)
#endif
  !output_flag=0 ! parameter
  !binary_flag=1
  call mpi_barrier(mpi_comm_world,ierr)
  do kk=1,2
     do isolve=1,num_solves
        ! flush cache
        do ii=1,dummy_size
           dummy(ii) = 1.235d3
        end do
#ifdef HAVE_PETSC
        if (kk==2) call PetscLogEventBegin(events(1),ierr)
#endif
        select case (problemType)
        case(0)
           call solve(grids,Apply_const_Lap,Apply_const_Lap,GSRB_const_Lap,res0,errors,nViters)
        case(1)
           if(mype==0) write(6,*) 'problemType 1 not implemented -- todo'
        end select
#ifdef HAVE_PETSC
        if (kk==2) call PetscLogEventEnd(events(1),ierr)
#endif
        if (kk==1) exit ! just do one solve to page everything in
        ! output run data and convergance data
        if(mype==0) then
           coarsest_grid = 0
           if (verbose .gt.1) write(6,*)'output: g(0).istop=',is_top(grids(0)),&
                'nfcycles=',nfcycles
           if (nfcycles .ne. 0 .and. .not. is_top(grids(0)) ) then
              write(iconv,900) isolve,-1.d0,-1.d0,errors(0)%uerror,&
                   errors(0)%graduerr,errors(0)%resid
              do ii=1,dom_max_grids-1
                 ! go from coarse to fine.  Need to pass through top empty 'coarse' grids
                 rateu = errors(ii-1)%uerror/errors(ii)%uerror
                 rategrad = errors(ii-1)%graduerr/errors(ii)%graduerr
                 if (verbose .gt.1) write(6,'(I5,A8,I2,A20,E14.6,A10,E14.6)') isolve,&
                      ': level ',ii, ': converg order u=', &
                      log(rateu)*log2r,', grad(u)=',log(rategrad)*log2r
                 write(iconv,900) isolve,log(rateu)*log2r,log(rategrad)*log2r,errors(ii)%uerror,&
                      errors(ii)%graduerr,errors(ii)%resid
                 
                 if( is_top(grids(ii)) ) then 
                    coarsest_grid = ii
                    exit ! grids and errors are in reverse order - this is just a counter
                 end if
              end do
           end if
 
           do ii=coarsest_grid,coarsest_grid+nViters-1 ! i is zero bases for 'errors'
              if (verbose.gt.2) write(6,'(A,I2,A,E14.6,A,E14.6)') 'driver: level',&
                   ii, ': error=',errors(ii)%uerror,&
                   ': resid=',errors(ii)%resid
              write(iconv,901) 'V:',isolve,0.d0,0.d0,errors(ii)%uerror,errors(ii)%graduerr,&
                   errors(ii)%resid
           end do
           
           ! general output
           ii = nfcycles*(coarsest_grid+1) + nViters - 1
           if (verbose.gt.1) write(6,888) 'solve ',isolve,' done |r|/|f|=',errors(ii)%resid/res0
           write(irun,700) isolve,errors(ii)%uerror,errors(ii)%graduerr,errors(ii)%resid/res0
        endif
     enddo
  end do
700 format (I7,6x,E14.6,E14.6,3x,E14.6)
900 format (I7,5x,E14.6,E14.6,3x,E14.6,E14.6,E14.6,E14.6)
901 format (A,I5,5x,E14.6,E14.6,3x,E14.6,E14.6,E14.6,E14.6)
888 format (A,I5,A,E14.6)
!!$  if(output_flag.eq.2) then
!!$     if(binary_flag.nq.0) then
!!$        if(mype==0) then
!!$           write(6,*) 'Binary output at',timeStep, ttot
!!$        endif
!!$        Call WriteSiloBinaryFileParallel(ux,timeStep)
!!$     endif
!!$  endif
#ifdef HAVE_PETSC
  call PetscFinalize(ierr)
#endif

  return
end subroutine driver
!-----------------------------------------------------------------
! get command line args, set problemType and fine grid sizes
!-----------------------------------------------------------------       
subroutine get_params(nprocs, nprocx, nprocy, nprocz, &
     nx, ny, nz, nxlocal, nylocal, nzlocal)
  use domain
  use mpistuff
  !use tags, only:msg_tag_inc
#ifndef HAVE_COMM_ARGS
  use f2kcli        ! command line args class
#endif
  implicit none
  integer,intent(out) :: nprocx,nprocy,nprocz,nx,ny,nz,nxlocal,nylocal,nzlocal
  integer,intent(in) :: nprocs
  !       
  CHARACTER(LEN=256) :: LINE
  CHARACTER(LEN=32)  :: CMD,CMD2
  INTEGER            :: NARG,IARG,fact,nsmooths
  integer*8 :: ti
  integer, parameter :: nargs=19
  double precision t2,t1,iarray(nargs)
  !       default input args: ?local, n?, ?procs
  nprocx=-1
  nprocy=-1
  nprocz=-1
  nx=-1
  ny=-1
  nz=-1
  nxlocal=-1
  nylocal=-1
  nzlocal=-1
  !       Problem type 
  !       0: (x^4-x^2)^D - Diri BCs
  !       1: bubble      - periodic domain
  problemType = 0
  ! solver parameters
  nvcycles = 0 ! pure FMG
  nfcycles = 1 ! pure FMG
  rtol = 1.d-6    ! only used for V-cycles
  ncycles = 1     ! v-cycles or w-cycles
  nfmgvcycles = 1 ! no interface for this (always 1)
  nsmooths = 2
  verbose = 1
  bot_min_size = 2 ! min size for bottom solver (grad get messed up with 1)
  mg_min_size = 8  ! or 16 ...
  num_solves = 1
  !msg_tag_inc = 0 ! better place to initilize this?
  ! read in args
  if( mype==0 ) then
     NARG = COMMAND_ARGUMENT_COUNT()
     DO IARG = 1,NARG,2
        CALL GET_COMMAND_ARGUMENT(IARG,CMD)
        CALL GET_COMMAND_ARGUMENT(IARG+1,CMD2)
        if( CMD == '-nx' ) READ (CMD2, '(I10)') nx
        if( CMD == '-ny' ) READ (CMD2, '(I10)') ny
        if( CMD == '-nz' ) READ (CMD2, '(I10)') nz
        
        if( CMD == '-nxloc' ) READ (CMD2, '(I10)') nxlocal
        if( CMD == '-nyloc' ) READ (CMD2, '(I10)') nylocal
        if( CMD == '-nzloc' ) READ (CMD2, '(I10)') nzlocal
        
        if( CMD == '-nxpes' ) READ (CMD2, '(I10)') nprocx
        if( CMD == '-nypes' ) READ (CMD2, '(I10)') nprocy
        if( CMD == '-nzpes' ) READ (CMD2, '(I10)') nprocz

        if( CMD == '-bot_min_size' ) READ (CMD2, '(I10)') bot_min_size
        if( CMD == '-mg_min_size' ) READ (CMD2, '(I10)') mg_min_size

        if( CMD == '-problem' ) READ (CMD2, '(I10)') problemType

        if( CMD == '-nvcycles' ) READ (CMD2, '(I10)') nvcycles
        if( CMD == '-nsmooths' ) READ (CMD2, '(I10)') nsmooths
        if( CMD == '-nfcycles' ) READ (CMD2, '(I10)') nfcycles
        if( CMD == '-rtol' )    READ (CMD2, '(E14.6)') rtol
        if( CMD == '-ncycle' )  READ (CMD2, '(I10)') ncycles

        if( CMD == '-verbose' )  READ (CMD2, '(I10)') verbose
        if( CMD == '-num_solves' )  READ (CMD2, '(I10)') num_solves
     END DO
#ifdef TWO_D
     nz = 1; nzlocal = 1; nprocz = 1;
#endif
     iarray(1) = nx
     iarray(2) = ny
     iarray(3) = nz
     
     iarray(4) = nxlocal
     iarray(5) = nylocal
     iarray(6) = nzlocal
     
     iarray(7) = nprocx
     iarray(8) = nprocy
     iarray(9) = nprocz
     
     iarray(10) = problemType
  
     iarray(11) = nvcycles
     iarray(12) = nfcycles
     iarray(13) = rtol
     iarray(14) = ncycles
     iarray(15) = nsmooths

     iarray(16) = verbose

     iarray(17) = bot_min_size
     iarray(18) = mg_min_size
     
     iarray(19) = num_solves

     call MPI_Bcast(iarray,nargs,MPI_DOUBLE_PRECISION, 0, &
          mpi_comm_world, ierr)
  else
     call MPI_Bcast(iarray,nargs,MPI_DOUBLE_PRECISION, 0, &
          mpi_comm_world, ierr)
     nx = int(iarray(1))
     ny = int(iarray(2))
     nz = int(iarray(3))
     
     nxlocal = int(iarray(4))
     nylocal = int(iarray(5))
     nzlocal = int(iarray(6))
     
     nprocx = int(iarray(7))
     nprocy = int(iarray(8))
     nprocz = int(iarray(9))
     
     problemType = int(iarray(10))

     nvcycles =    int(iarray(11))
     nfcycles =    int(iarray(12))
     rtol     =        iarray(13) 
     ncycles  =    int(iarray(14))
     nsmooths =    int(iarray(15))

     verbose  =    int(iarray(16))

     bot_min_size  = int(iarray(17))
     mg_min_size   = int(iarray(18))

     num_solves = int(iarray(19))
  endif
  if (nfcycles .ne. 0) nfcycles = 1 ! only valid 0,1
  if (nvcycles .lt. 0) nvcycles = 0 ! only valid 0,1
  if (nfcycles+nvcycles .le. 0) print *, 'no solver iterations!!!'
  nsmoothsup = nsmooths
  nsmoothsdown = nsmooths
  ncoarsesolveits = bot_min_size*bot_min_size

  ! get ijk procs
  t1=nprocs
#ifdef TWO_D
  nprocz = 1
  nz = 1
  nzlocal = 1
  if (nprocx == -1) then
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
     nprocx=fact*floor(dsqrt(t1/fact))
     nprocy=t1/nprocx
     if (nprocs .ne. nprocz*nprocx*nprocy) then
        write(6,*) 'INCORRECT NP:', NPROCS,'nprocx=',nprocx,&
             'nprocy=',nprocy,'nprocz=',nprocz,'nprocs=',nprocs
        stop
     endif
  else
     nprocy=t1/nprocx
     if(nprocy<1)stop'too many nprocx???'
     fact = nprocx/nprocy
     if( fact < 1 ) then
        write(6,*) 'nprocx must be greater than nprocy', NPROCS,&
       'nprocx=',nprocx,&
       'nprocy=',nprocy,'nprocz=',nprocz,'nprocs=',nprocs
        stop
     endif
  endif
#else
  if (nprocx == -1) then
     t1 = nprocs
     if(nprocs==1) then
        nprocx = 1; nprocy = 1; nprocz = 1;
     elseif(nx==ny .and. ny==nz)then
        nprocx = floor(t1**.3333333334d0)
        nprocy = nprocx
        nprocz = nprocs/(nprocx*nprocy)
     elseif(nz==ny)then
        t1 = ny*nprocs; t2 = nx; t1 = t1/t2 
        nprocy = floor(t1**.333333334d0)
        nprocz = nprocy
        nprocx = nprocs/(nprocy*nprocz) 
     elseif(nz==nx)then
        t1 = nz*nprocs; t2 = ny; t1 = t1/t2 
        nprocx = floor(t1**.333333334d0)
        nprocz = nprocx
        nprocy = nprocs/(nprocx*nprocz) 
     elseif(nx==ny)then
        t1 = ny*nprocs; t2 = nz; t1 = t1/t2 
        nprocy = floor(t1**.333333334d0)
        nprocx = nprocy
        nprocz = nprocs/(nprocy*nprocx) 
     else
        stop 'need to specify x/y/z proc'
     endif
  endif
#endif     
  if (nprocs .ne. nprocz*nprocx*nprocy) then
     write(6,*) 'INCORRECT NP: nprocx=',nprocx, &
          'nprocy=',nprocy,'nprocz=',nprocz,'nprocs=',nprocs
     stop
  endif
  ! x
  if (nxlocal.eq.-1 .and. nx.eq.-1) nxlocal = 64        ! default
  if (nxlocal.eq.-1 .and. nx.ne.-1) nxlocal = nx/nprocx 
  if (nxlocal.ne.-1 .and. nx.eq.-1) nx = nxlocal*nprocx 
  if (nx .ne. nxlocal*nprocx) stop 'nx .ne. nxlocal*nprocx'
  ! y
  if (nylocal.eq.-1 .and. ny.eq.-1) nylocal = nxlocal   ! default square
  if (nylocal.eq.-1 .and. ny.ne.-1) nylocal = ny/nprocy 
  if (nylocal.ne.-1 .and. ny.eq.-1) ny = nylocal*nprocy 
  if (ny .ne. nylocal*nprocy) stop 'ny .ne. nylocal*nprocy'
  ! z
  if (nzlocal.eq.-1 .and. nz.eq.-1) nzlocal = nxlocal   ! default square
  if (nzlocal.eq.-1 .and. nz.ne.-1) nzlocal = nz/nprocz 
  if (nzlocal.ne.-1 .and. nz.eq.-1) nz = nzlocal*nprocz 
  if (nz .ne. nzlocal*nprocz) stop 'nz .ne. nzlocal*nprocz'
  ! print topology
  if(mype==0 .and. verbose .gt. 3) then
     write(6,*) '[',mype,'] nx =',nx,     ',ny= ',ny     ,',nz= ',nz
     write(6,*) '[',mype,'] npx=',nprocx, ',npy=',nprocy, ',npz=',nprocz
     write(6,*) '[',mype,'] nxl=',nxlocal,',nyl=',nylocal,',nzl=',nzlocal
  endif
end subroutine get_params
!-----------------------------------------------------------------------
subroutine formExactU(exact,g)
  use GridModule
  use mpistuff
  use domain
  implicit none
  type(proc_patch),intent(in):: g
  double precision,intent(out)::exact(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar) 
  !
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(7),ierr)
#endif
  ig = g%iglobalx
  jg = g%iglobaly
  kg = g%iglobalz
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
              !write (6,'(A18,E14.6,E14.6,E14.6)') "formExactU: coord=",coord
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
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(7),ierr)
#endif
  call SetBCs(exact,g)
end subroutine formExactU
!-----------------------------------------------------------------------
subroutine formGradError(exactu,ux,g,gerr,order)
  use GridModule
  use domain
  use mpistuff
  implicit none
  type(proc_patch),intent(in):: g
  double precision,intent(in)::exactu(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar) 
  double precision,intent(in)::ux    (g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar) 
  double precision,intent(out)::gerr
  integer,intent(in):: order
  !
  integer:: ii,jj,kk
  double precision::grad1(3),grad2(3)
  double precision::gradError(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  double precision norm
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(7),ierr)
#endif
  select case (problemType)
  case(0)
     do kk=1,g%kmax
        do jj=1,g%jmax
           do ii=1,g%imax
              grad1(1) = 0.5d0*(ux(ii+1,jj,kk,1)-ux(ii-1,jj,kk,1))/g%dxg
              grad1(2) = 0.5d0*(ux(ii,jj+1,kk,1)-ux(ii,jj-1,kk,1))/g%dyg
              grad2(1) = 0.5d0*(exactu(ii+1,jj,kk,1)-exactu(ii-1,jj,kk,1))/g%dxg
              grad2(2) = 0.5d0*(exactu(ii,jj+1,kk,1)-exactu(ii,jj-1,kk,1))/g%dyg
              
#ifndef TWO_D
              grad1(3) = 0.5d0*(ux(ii,jj,kk+1,1)-ux(ii,jj,kk-1,1))/g%dzg
              grad2(3) = 0.5d0*(exactu(ii,jj,kk+1,1)-exactu(ii,jj,kk-1,1))/g%dzg
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
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(7),ierr)
#endif
end subroutine formGradError
!-----------------------------------------------------------------------
subroutine FormRHS(rhs,g)
  use GridModule
  use domain
  use mpistuff
  implicit none
  type(proc_patch),intent(in):: g
  double precision,intent(out)::rhs(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(7),ierr)
#endif 
  ig = g%iglobalx
  jg = g%iglobaly
  kg = g%iglobalz
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
#ifndef TWO_D
              coord(3) = zl+(kg+kk-1)*g%dzg-0.5*g%dzg          
#endif
              !write (6,'(A15,E14.6,E14.6,E14.6)') "FormRHS: coord=",coord
              x2 = coord*coord
              a = x2*(x2-1.);              
              b = 12.*x2-2.;
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
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(7),ierr)
#endif
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
  use error_data_module, only:error_norm
  implicit none
  if(mype==0) then
     open(irun,file='Run.history', form='formatted')
      write(irun,  '(A,I1,A)') 'solve # error_',error_norm,': u           grad(u)          |f-Au|_2/|f|_2'
     write(irun,*) '-----------------------------------------------------------'
     open(iconv,file='Convergence.history', form='formatted')
      write(iconv, '(A,I1,A)')'solve # order (',error_norm, ' norm): u   grad(u)          error: u      grad(u)       |f-Au|_2'
     write(iconv,*)'-------------------------------------------------------------------------------'
  endif
  return
end subroutine InitIO
!-----------------------------------------------------------------------
subroutine FinalizeIO()
  use iounits
  use mpistuff, only:mype
  implicit none
  ! Close all open files
  if(mype==0) then
     close(irun)
     close(iconv)
  endif
  return
end subroutine FinalizeIO
