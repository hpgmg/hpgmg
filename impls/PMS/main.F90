!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
program main
  use grid_module
  use domain, only:periodic
  use pms, only:verbose,dom_max_grids,dom_max_sr_grids,problemType
  use mpistuff
  use tags
  !use discretization
  implicit none
  integer:: rank,comm3d,ndims
  integer:: NPeAxis(3),iPeAxis(3)
  type(ipoint):: nloc,npe
  type(crs_patcht):: gc(0:dom_max_grids-1)
  type(sr_patcht):: gsr(-dom_max_sr_grids:0)

  ! global constructor stuff
  MSG_RESTRICT_TAG = 90 ! this incs slower than the XCH
  MSG_XCH_XLOW_TAG=100
  MSG_XCH_XHI_TAG=200
  MSG_XCH_YLOW_TAG=300
  MSG_XCH_YHI_TAG=400
  MSG_XCH_ZLOW_TAG=500
  MSG_XCH_ZHI_TAG=600

  ! MPI init - globals
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, mype, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, mpisize, ierr)

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

  ! get fine grid parameters, problem type, and global variables
  call get_params(npe,nloc)

  ! Initialize the fine grid communicator, set periodicity
  ndims = 3
  NPeAxis(1) = npe%i
  NPeAxis(2) = npe%j
  NPeAxis(3) = npe%k
  periodic(1) = .false.
  periodic(2) = .false.
  periodic(3) = .false.
#ifdef XPERIODIC
  periodic(1) = .true.
#endif  
#ifdef YPERIODIC
  periodic(2) = .true.
#endif  
#ifdef ZPERIODIC
  periodic(3) = .true.
#endif
  call MPI_cart_create(MPI_COMM_WORLD,ndims,NPeAxis,periodic,&
       .true.,comm3d,ierr) ! reoder is ignored???
  call MPI_comm_rank(comm3d,rank,ierr) ! new rank?
  if ((mype==0.and.verbose.gt.3).or.mype.ne.rank) then
     if (verbose.gt.0) write(6,*) 'main: warning old mype=',mype,', comm3D new rank = ',rank
  endif
  mype=rank

  call MPI_Cart_Coords(comm3d,rank,ndims,iPeAxis,ierr)
  iPeAxis(1) = iPeAxis(1) + 1 ! one based
  iPeAxis(2) = iPeAxis(2) + 1
  iPeAxis(3) = iPeAxis(3) + 1
  
  ! run it
  call SetupDomain(problemType) ! sets size of domain

  ! creates grids meta data, comms, sets global topo vars
  call new_grids(gc,gsr,NPeAxis,iPeAxis,nloc,comm3d)

  ! do it 
  call driver(gc,gsr)
  
  ! end it
  call destroy_grids_private(gc,gsr)
  call FinalizeIO()
  call MPI_Finalize(ierr)
end program MAIN
!-----------------------------------------------------------------
!  driver - the PETSc main
!-----------------------------------------------------------------
subroutine driver(gc,gsr)
  use iounits
  use mpistuff 
  use pms
  use grid_module
  use error_data_module
  implicit none
  !
  type(crs_patcht),intent(out):: gc(0:dom_max_grids-1)
  type(sr_patcht),intent(out):: gsr(-dom_max_sr_grids:0)
  !
  integer:: isolve,ii,nViters,kk
  integer,parameter:: cache_flusher_size=1 !1024*1024 ! cache flush array size
  double precision:: cache_flusher(cache_flusher_size),res0,rateu,rategrad
  type(error_data)::errors(0:ncgrids-1+nvcycles+dom_max_sr_grids)
  !external Apply_const_Lap,GS_RB_const_Lap,Jacobi_const_Lap
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
  call mpi_barrier(mpi_comm_world,ierr)
  do kk=1,1
     do isolve=1,num_solves
        ! flush cache
        do ii=1,cache_flusher_size
           cache_flusher(ii) = 1.235d3
        end do
#ifdef HAVE_PETSC
        if (kk==1) call PetscLogEventBegin(events(1),ierr)
#endif
        select case (problemType)
        case(0)
           call solve(gc,gsr,Apply_const_Lap,Apply_const_Lap,GS_RB_const_Lap,res0,errors,nViters)
        case(1)
           if (mype==0) write(6,*) 'problemType 1 not implemented -- todo'
        end select
#ifdef HAVE_PETSC
        if (kk==1) call PetscLogEventEnd(events(1),ierr)
#endif
        !if (kk==1) exit ! just do one solve to page everything in
        ! output run data and convergance data
        if (mype==0) then
           if (nfcycles.ne.0 .and. .not.is_top(gc(0)) ) then
              write(iconv,900) isolve,-1.d0,-1.d0,errors(0)%uerror,&
                   -1.0,errors(0)%resid
              do ii=1,dom_max_grids-1+nsr
                 ! go from coarse to fine.  Need to pass through top empty 'coarse' grids
                 rateu = errors(ii-1)%uerror/errors(ii)%uerror
                 rategrad = 1.0
                 if (verbose .gt.0) write(6,'(A,I2,A,E14.6)') &
                      'driver: level ',ii, ': converg order u=', &
                      log(rateu)*log2r
                 write(iconv,900) isolve,log(rateu)*log2r,log(rategrad)*log2r,errors(ii)%uerror,&
                      -1.0,errors(ii)%resid
                 if (ii.gt.nsr.and.is_top(gc(ii-nsr))) exit 
              end do
           end if
           ! print V cycles
           if (nViters>0) then
              do ii=0,nViters-1 ! i is zero bases for 'errors'
                 if (verbose.gt.0) write(6,'(A,I2,A,E14.6,A,E14.6)') 'driver: V cycle',&
                      ii, ': error=',errors(ii)%uerror,&
                      ': resid=',errors(ii)%resid
                 write(iconv,901) 'V:',isolve,0.d0,0.d0,errors(ii)%uerror,-1.0,&
                      errors(ii)%resid
              end do
           end if
           ! general output
           ii = ncgrids + nViters - 1 + nsr
           if (verbose.gt.1) write(6,888) 'solve ',isolve,' done |r|/|f|=',errors(ii)%resid/res0
           write(irun,700) isolve,errors(ii)%uerror,-1.0,errors(ii)%resid/res0
        endif
     enddo
  end do
700 format (I7,6x,E14.6,E14.6,3x,E14.6)
900 format (I7,5x,E14.6,E14.6,3x,E14.6,E14.6,E14.6,E14.6)
901 format (A,I5,5x,E14.6,E14.6,3x,E14.6,E14.6,E14.6,E14.6)
888 format (A,I5,A,E14.6)
#ifdef HAVE_PETSC
  call PetscFinalize(ierr)
#endif
  return
end subroutine driver
!-----------------------------------------------------------------
! get command line args, set problemType and fine grid sizes
!-----------------------------------------------------------------       
subroutine get_params(npe,nloc)
  use pms
  use mpistuff
  use base_data_module
  use error_data_module,only:error_norm
#ifndef HAVE_COMM_ARGS
  use f2kcli        ! command line args class
#endif

  implicit none
  type(ipoint),intent(out)::npe,nloc
  !
  type(ipoint)::nglob
  CHARACTER(LEN=256) :: LINE
  CHARACTER(LEN=32)  :: CMD,CMD2
  INTEGER            :: NARG,IARG,fact,nsmooths
  integer*8 :: ti
  integer, parameter :: nargs=26
  double precision t2,t1,iarray(nargs)
  ! default input args: ?local, n?, ?pes
  npe%i=-1
  npe%j=-1
  npe%k=-1
  nglob%i=-1
  nglob%j=-1
  nglob%k=-1
  nloc%i=-1
  nloc%j=-1
  nloc%k=-1
  !       Problem type 
  !       0: (x^4-x^2)^D - Diri BCs
  !       1: bubble      - periodic domain???
  problemType = 0
  ! solver parameters
  error_norm = 3 ! inf==3
  nvcycles = 0 ! pure FMG
  nfcycles = 1 ! pure FMG
  rtol = 1.d-5    ! only used for V-cycles
  ncycles = 1     ! v-cycles or w-cycles
  nfmgvcycles = 1 ! no interface for this (always 1)
  nsmooths = 1
  nsmoothsup = -1
  nsmoothsdown = -1
  nsmoothsfmg = -1
  verbose = 1
  bot_min_sz = 2 ! min size for bottom solver (grad get messed up with 1)
  pe_min_sz = 8  ! or 16 ...
  num_solves = 1
  sr_max_loc_sz = 0 ! no SR by default
  sr_base_bufsz = 4
  sr_bufsz_inc = 0
  !msg_tag_inc = 0 ! better place to initilize this?
  ! read in args
  if ( mype==0 ) then
     NARG = COMMAND_ARGUMENT_COUNT()
     DO IARG = 1,NARG,2
        CALL GET_COMMAND_ARGUMENT(IARG,CMD)
        CALL GET_COMMAND_ARGUMENT(IARG+1,CMD2)
        if ( CMD == '-nx' ) READ (CMD2, '(I10)') nglob%i
        if ( CMD == '-ny' ) READ (CMD2, '(I10)') nglob%j
        if ( CMD == '-nz' ) READ (CMD2, '(I10)') nglob%k
        
        if ( CMD == '-nxloc' ) READ (CMD2, '(I10)') nloc%i
        if ( CMD == '-nyloc' ) READ (CMD2, '(I10)') nloc%j
        if ( CMD == '-nzloc' ) READ (CMD2, '(I10)') nloc%k
        
        if ( CMD == '-nxpe' ) READ (CMD2, '(I10)') npe%i
        if ( CMD == '-nype' ) READ (CMD2, '(I10)') npe%j
        if ( CMD == '-nzpe' ) READ (CMD2, '(I10)') npe%k

        if ( CMD == '-bot_min_sz' ) READ (CMD2, '(I10)') bot_min_sz
        if ( CMD == 'pe_min_sz' ) READ (CMD2, '(I10)') pe_min_sz

        if ( CMD == '-problem' ) READ (CMD2, '(I10)') problemType

        if ( CMD == '-nvcycles' ) READ (CMD2, '(I10)') nvcycles
        if ( CMD == '-nsmooths' ) READ (CMD2, '(I10)') nsmooths
        if ( CMD == '-nfcycles' ) READ (CMD2, '(I10)') nfcycles
        if ( CMD == '-rtol' )    READ (CMD2, '(E14.6)') rtol
        if ( CMD == '-ncycle' )  READ (CMD2, '(I10)') ncycles

        if ( CMD == '-verbose' )  READ (CMD2, '(I10)') verbose
        if ( CMD == '-num_solves' )  READ (CMD2, '(I10)') num_solves

        if ( CMD == '-sr_max_loc_sz' )  READ (CMD2, '(I10)') sr_max_loc_sz

        if ( CMD == '-sr_base_bufsz' )  READ (CMD2, '(I10)') sr_base_bufsz
        if ( CMD == '-sr_bufsz_inc' )  READ (CMD2, '(I10)') sr_bufsz_inc

        if ( CMD == '-nsmoothsup' ) READ (CMD2, '(I10)') nsmoothsup
        if ( CMD == '-nsmoothsdown' ) READ (CMD2, '(I10)') nsmoothsdown
        if ( CMD == '-nsmoothsfmg' ) READ (CMD2, '(I10)') nsmoothsfmg

        if ( CMD == '-error_norm' ) READ (CMD2, '(I10)') error_norm
     END DO
#ifdef TWO_D
     nglob%k = 1; nloc%k = 1; npe%k = 1;
#endif
     iarray(1) = nglob%i
     iarray(2) = nglob%j
     iarray(3) = nglob%k
     
     iarray(4) = nloc%i
     iarray(5) = nloc%j
     iarray(6) = nloc%k
     
     iarray(7) = npe%i
     iarray(8) = npe%j
     iarray(9) = npe%k
     
     iarray(10) = problemType
  
     iarray(11) = nvcycles
     iarray(12) = nfcycles
     iarray(13) = rtol
     iarray(14) = ncycles
     iarray(15) = nsmooths

     iarray(16) = verbose

     iarray(17) = bot_min_sz
     iarray(18) = pe_min_sz
     
     iarray(19) = num_solves

     iarray(20) = sr_max_loc_sz
     iarray(21) = sr_base_bufsz
     iarray(22) = sr_bufsz_inc

     iarray(23) = nsmoothsup
     iarray(24) = nsmoothsdown
     iarray(25) = nsmoothsfmg

     iarray(26) = error_norm

     call MPI_Bcast(iarray,nargs,MPI_DOUBLE_PRECISION, 0, &
          mpi_comm_world, ierr)
  else
     call MPI_Bcast(iarray,nargs,MPI_DOUBLE_PRECISION, 0, &
          mpi_comm_world, ierr)
     nglob%i = int(iarray(1))
     nglob%j = int(iarray(2))
     nglob%k = int(iarray(3))
     
     nloc%i = int(iarray(4))
     nloc%j = int(iarray(5))
     nloc%k = int(iarray(6))
     
     npe%i = int(iarray(7))
     npe%j = int(iarray(8))
     npe%k = int(iarray(9))
     
     problemType = int(iarray(10))

     nvcycles =    int(iarray(11))
     nfcycles =    int(iarray(12))
     rtol     =        iarray(13) 
     ncycles  =    int(iarray(14))
     nsmooths =    int(iarray(15))

     verbose  =    int(iarray(16))

     bot_min_sz  = int(iarray(17))
     pe_min_sz   = int(iarray(18))

     num_solves = int(iarray(19))
     
     sr_max_loc_sz = int(iarray(20))
     sr_base_bufsz = int(iarray(21))
     sr_bufsz_inc = int(iarray(22))
     
     nsmoothsup = int(iarray(23))
     nsmoothsdown = int(iarray(24))
     nsmoothsfmg = int(iarray(25))

     error_norm = int(iarray(26))
  endif
  if (nfcycles .ne. 0) nfcycles = 1 ! only valid 0,1
  if (nvcycles .lt. 0) nvcycles = 0 ! only valid 0,1
  if (nfcycles+nvcycles .le. 0) print *, 'no solver iterations!!!'
  ! general constructor stuff
  if (nsmoothsup==-1) nsmoothsup = nsmooths
  if (nsmoothsdown==-1) nsmoothsdown = nsmooths
  if (nsmoothsfmg==-1) nsmoothsfmg = 0 ! default
  ncoarsesolveits = bot_min_sz*bot_min_sz
  
  if (verbose.gt.3) then
     print *,'problemType=',problemType
     print *,'nvcycles=',nvcycles
     print *,'nfcycles=',nfcycles
     print *,'rtol=',rtol
     print *,'ncycles=',ncycles
     print *,'bot_min_sz=',bot_min_sz
     print *,'pe_min_sz=',pe_min_sz
     print *,'num_solves=',num_solves
     print *,'sr_max_loc_sz=',sr_max_loc_sz
     print *,'sr_base_bufsz=',sr_base_bufsz
     print *,'sr_bufsz_inc=',sr_bufsz_inc
     print *,'nsmoothsup=',nsmoothsup
     print *,'nsmoothsdown=',nsmoothsdown
     print *,'nsmoothsfmg=',nsmoothsfmg
     print *,'error_norm=',error_norm
  end if

  ! get ijk pes
  t1 = mpisize
#ifdef TWO_D
  npe%k = 1
  nglob%k = 1
  nloc%k = 1
  if (npe%i == -1) then
     if ( abs(dsqrt(t1) - floor(dsqrt(t1))) > 0.d0 ) then
        if (abs(dsqrt(t1/2.d0)-floor(dsqrt(t1/2.d0)))>0.d0)then
           if (abs(dsqrt(t1/4.d0)-floor(dsqrt(t1/4.d0)))>0.d0)then
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
     npe%i=fact*floor(dsqrt(t1/fact))
     npe%j=t1/npe%i
     if (mpisize .ne. npe%k*npe%i*npe%j) then
        write(6,*) 'INCORRECT NP:', NPES,'npe%i=',npe%i,&
             'npe%j=',npe%j,'npe%k=',npe%k,'npes=',mpisize
        stop
     endif
  else (npe%j == -1) then
     npe%j=t1/npe%i
     if (npe%j<1)stop'too many npe%i???'
     fact = npe%i/npe%j
     if ( fact < 1 ) then
        write(6,*) 'npe%i must be greater than npe%j', NPES,&
       'npe%i=',npe%i,&
       'npe%j=',npe%j,'npe%k=',npe%k,'npes=',mpisize
        stop
     endif
  endif
#else
  if (npe%i == -1) then
     t1 = mpisize
     if (mpisize==1) then ! one pe
        npe%i = 1; npe%j = 1; npe%k = 1;
     elseif (nglob%i==nglob%j .and. nglob%j==nglob%k)then ! square
        npe%i = floor(t1**.3333333334d0)
        npe%j = npe%i
        npe%k = mpisize/(npe%i*npe%j)
     elseif (nglob%k==nglob%j)then ! z==y
        t1 = nglob%j*mpisize; t2 = nglob%i; t1 = t1/t2 
        npe%j = floor(t1**.333333334d0)
        npe%k = npe%j
        npe%i = mpisize/(npe%j*npe%k) 
     elseif (nglob%k==nglob%i)then
        t1 = nglob%k*mpisize; t2 = nglob%j; t1 = t1/t2 
        npe%i = floor(t1**.333333334d0)
        npe%k = npe%i
        npe%j = mpisize/(npe%i*npe%k) 
     elseif (nglob%i==nglob%j)then
        t1 = nglob%j*mpisize; t2 = nglob%k; t1 = t1/t2 
        npe%j = floor(t1**.333333334d0)
        npe%i = npe%j
        npe%k = mpisize/(npe%j*npe%i) 
     else
        stop 'need to specify x/y/z pe'
     endif
  endif
#endif     
  if (mpisize .ne. npe%k*npe%i*npe%j) then
     write(6,*) 'INCORRECT NP: npe%i=',npe%i, &
          'npe%j=',npe%j,'npe%k=',npe%k,'npes=',mpisize
     stop
  endif
  ! x
  if (nloc%i.eq.-1 .and. nglob%i.eq.-1) nloc%i = 64        ! default
  if (nloc%i.eq.-1 .and. nglob%i.ne.-1) nloc%i = nglob%i/npe%i 
  if (nloc%i.ne.-1 .and. nglob%i.eq.-1) nglob%i = nloc%i*npe%i 
  if (nglob%i .ne. nloc%i*npe%i) stop 'nx .ne. nloc%i*npe%i'
  ! y
  if (nloc%j.eq.-1 .and. nglob%j.eq.-1) nloc%j = nloc%i   ! default square
  if (nloc%j.eq.-1 .and. nglob%j.ne.-1) nloc%j = nglob%j/npe%j 
  if (nloc%j.ne.-1 .and. nglob%j.eq.-1) nglob%j = nloc%j*npe%j 
  if (nglob%j .ne. nloc%j*npe%j) stop 'ny .ne. nylocal*npe%j'
  ! z
  if (nloc%k.eq.-1 .and. nglob%k.eq.-1) nloc%k = nloc%i   ! default square
  if (nloc%k.eq.-1 .and. nglob%k.ne.-1) nloc%k = nglob%k/npe%k 
  if (nloc%k.ne.-1 .and. nglob%k.eq.-1) nglob%k = nloc%k*npe%k 
  if (nglob%k .ne. nloc%k*npe%k) stop 'nz .ne. nloc%k*npe%k'
  ! print topology
  if (mype==0 .and. verbose .gt. 3) then
     write(6,*) '[',mype,'] nx =',nglob%i,',ny= ',nglob%j,',nz= ',nglob%k
     write(6,*) '[',mype,'] npx=',npe%i,  ',npy=',npe%j,  ',npz=',npe%k
     write(6,*) '[',mype,'] nxl=',nloc%i, ',nyl=',nloc%j, ',nzl=',nloc%k
  endif
end subroutine get_params

!-----------------------------------------------------------------
subroutine SetupDomain(problemType)
  use domain
  implicit none
  integer::problemType
  !
  !     Problem type 
  !     0: (x^4-x^2)
  !     1: bubble
  !
  if (problemType.eq.0) then
     xl=0.D0
     xr=1.D0
     yl=0.D0
     yr=1.D0
     zl=0.D0
     zr=1.D0
  else if (problemType.eq.1) then
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
  use error_data_module,only:error_norm
  implicit none
  if (mype==0) then
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
subroutine finalizeio()
  use iounits
  use mpistuff, only:mype
  implicit none
  ! Close all open files
  if (mype==0) then
     close(irun)
     close(iconv)
  endif
  return
end subroutine finalizeio
!-----------------------------------------------------------------
subroutine destroy_grids_private(cg,gsr)
  use grid_module
  use mpistuff
  use pms, only:pe_min_sz,dom_max_sr_grids,dom_max_grids,nsr,ncgrids
  implicit none
  type(crs_patcht)::cg(0:dom_max_grids-1)
  type(sr_patcht)::gsr(-dom_max_sr_grids:0)
 
  integer :: n
  
  do n=1,ncgrids
     if (&
          (cg(n-1)%p%max%hi%i) > pe_min_sz .and. &
          (cg(n-1)%p%max%hi%j) > pe_min_sz &
#ifdef TWO_D
          ) then
#else
        .and. (cg(n-1)%p%max%hi%k) > pe_min_sz ) then
#endif
        ! normal reduction, no split yet -- nothing to destroy
     else 
        ! split & not same comms as last
        if (cg(n)%t%npe%i .ne. cg(n-1)%t%npe%i ) then
           call MPI_COMM_free(cg(n)%t%comm,ierr)
           call MPI_COMM_free(cg(n)%t%loc_comm,ierr) 
           call MPI_comm_free(cg(n)%t%comm3d, ierr)
        else
           exit
        end if
     end if
  enddo
  return
end subroutine destroy_grids_private
!-----------------------------------------------------------------
subroutine new_grids_private(cg,gsr,NPeAxis,iPeAxis,nloc,comm3d)
  use grid_module
  use mpistuff
  use pms
  use domain
  use discretization,only:nsg
  implicit none
  type(ipoint)::nloc
  integer,intent(in)::comm3d
  integer,intent(in)::NPeAxis(3),iPeAxis(3)
  type(crs_patcht),intent(out)::cg(0:dom_max_grids-1)
  type(sr_patcht),intent(out)::gsr(-dom_max_sr_grids:0)
  ! 
  integer::n,ndims,ii,jj,kk,srbf,pe_id,rank,npes(3)

  ndims = 3
  !----------------------------------------------------------------------
  ! set up finest (coarse) grid
  cg(0)%p%dx%i=(xr-xl)/(nloc%i*NPeAxis(1))
  cg(0)%p%dx%j=(yr-yl)/(nloc%j*NPeAxis(2))
  cg(0)%p%dx%k=(zr-zl)/(nloc%k*NPeAxis(3))

  cg(0)%t%comm3d = comm3d
  cg(0)%t%comm = comm3d
  cg(0)%t%loc_comm = mpi_comm_null 
  !     left, ...
  call MPI_Cart_Shift(comm3D,0,1,cg(0)%t%left,cg(0)%t%right,ierr)
  call MPI_Cart_Shift(comm3D,1,1,cg(0)%t%bottom,cg(0)%t%top,ierr)
  call MPI_cart_Shift(comm3D,2,1,cg(0)%t%behind,cg(0)%t%forward,ierr)
  if (mype==0.and.verbose.gt.4) write(6,*)'[',mype,'] l,r,t,b=',&
       cg(0)%t%left,cg(0)%t%right,cg(0)%t%top,&
       cg(0)%t%bottom,cg(0)%t%behind,cg(0)%t%forward
  !       ipex ...	
  cg(0)%t%ipe%i = iPeAxis(1)
  cg(0)%t%ipe%j = iPeAxis(2)
  cg(0)%t%ipe%k = iPeAxis(3)
  
  cg(0)%t%npe%i = NPeAxis(1) 
  cg(0)%t%npe%j = NPeAxis(2) 
  cg(0)%t%npe%k = NPeAxis(3) 
  ! 
  cg(0)%p%max%hi%i=nloc%i
  cg(0)%p%max%hi%j=nloc%j
  cg(0)%p%max%lo%i=1
  cg(0)%p%max%lo%j=1
#ifdef TWO_D
  cg(0)%p%max%hi%k=1 
  cg(0)%p%max%lo%k=1
#else 
  cg(0)%p%max%hi%k=nloc%k
  cg(0)%p%max%lo%k=1
#endif
  if (mype==0.and.verbose.gt.3) then
     write(6,*) '[',mype,'] level ',0,', nxl=',cg(0)%p%max%hi%i,'npx=',cg(0)%t%npe%i
  end if
  !----------------------------------------------------------------------
  ! do coarse grids (non-sr). start at 1 as 0 was just done
  ncgrids = 1
  do n=1,dom_max_grids-1
     cg(n)%p%dx%i = cg(n-1)%p%dx%i*mg_ref_ratio
     cg(n)%p%dx%j = cg(n-1)%p%dx%j*mg_ref_ratio
#ifdef TWO_D
     cg(n)%p%dx%k = 0.d0
#else 
     cg(n)%p%dx%k = cg(n-1)%p%dx%k*mg_ref_ratio
#endif
     ! take care of reductions
     if ( is_top(cg(n-1)) ) then
        ! all done - clear rest of grids
        do ii=n,dom_max_grids-1
           cg(ii)%t%npe%i = 0 
           cg(ii)%t%npe%j = 0 
           cg(ii)%t%npe%k = 0 
           cg(ii)%p%max%hi%i = 0 ! zero size 
           cg(ii)%p%max%hi%j = 0 
           cg(ii)%p%max%hi%k = 0 
           cg(ii)%p%max%lo%i = 1 
           cg(ii)%p%max%lo%j = 1 
           cg(ii)%p%max%lo%k = 1 
           ! malloc of next grid in MGV, etc
           cg(ii)%p%all%lo%i=1 ! zero size 
           cg(ii)%p%all%hi%i=0
           cg(ii)%p%all%lo%j=1
           cg(ii)%p%all%hi%j=0
           cg(ii)%p%all%lo%k=1
           cg(ii)%p%all%hi%k=0
           
           cg(ii)%t%loc_comm = mpi_comm_null
           cg(ii)%t%comm3d = mpi_comm_null
           cg(ii)%t%comm = mpi_comm_null
        end do
        exit
     else if (&
          (cg(n-1)%p%max%hi%i)>pe_min_sz .and. &
          (cg(n-1)%p%max%hi%j)>pe_min_sz &
#ifndef TWO_D
          .and.&
          (cg(n-1)%p%max%hi%k)>pe_min_sz &
#endif
          ) then
        !     normal reduction, no split yet	      
        cg(n)%p%max%hi%i=cg(n-1)%p%max%hi%i/mg_ref_ratio
        cg(n)%p%max%hi%j=cg(n-1)%p%max%hi%j/mg_ref_ratio
        cg(n)%p%max%lo%i=1
        cg(n)%p%max%lo%j=1
        cg(n)%p%max%lo%k=1
#ifndef TWO_D
        cg(n)%p%max%hi%k=cg(n-1)%p%max%hi%k/mg_ref_ratio
#else
        cg(n)%p%max%hi%k=1
#endif
        cg(n)%t%comm3d = comm3d
        cg(n)%t%comm = MPI_COMM_WORLD
        cg(n)%t%loc_comm = mpi_comm_null
        cg(n)%t%ipe%i = cg(n-1)%t%ipe%i
        cg(n)%t%ipe%j = cg(n-1)%t%ipe%j
        cg(n)%t%ipe%k = cg(n-1)%t%ipe%k
        cg(n)%t%npe%i = cg(n-1)%t%npe%i
        cg(n)%t%npe%j = cg(n-1)%t%npe%j
        cg(n)%t%npe%k = cg(n-1)%t%npe%k
        
        cg(n)%t%left = cg(n-1)%t%left
        cg(n)%t%right = cg(n-1)%t%right
        cg(n)%t%top = cg(n-1)%t%top
        cg(n)%t%bottom = cg(n-1)%t%bottom
        cg(n)%t%forward = cg(n-1)%t%forward
        cg(n)%t%behind = cg(n-1)%t%behind
        
     else if (cg(n-1)%t%npe%k==1 .or. cg(n-1)%t%npe%j==1 &
#ifndef TWO_D
          .or. cg(n-1)%t%npe%i==1 ) then ! one pe reduction, copy comm stuff
#else
        ) then
#endif
        cg(n)%t%npe%i = cg(n-1)%t%npe%i
        cg(n)%t%npe%j = cg(n-1)%t%npe%j
        cg(n)%t%ipe%i = cg(n-1)%t%ipe%i
        cg(n)%t%ipe%j = cg(n-1)%t%ipe%j
        cg(n)%t%npe%k = cg(n-1)%t%npe%k
        cg(n)%t%ipe%k = cg(n-1)%t%ipe%k
        
        cg(n)%p%max%hi%i=cg(n-1)%p%max%hi%i/mg_ref_ratio
        cg(n)%p%max%hi%j=cg(n-1)%p%max%hi%j/mg_ref_ratio
        cg(n)%p%max%lo%i=1
        cg(n)%p%max%lo%j=1
        cg(n)%p%max%lo%k=1
#ifndef TWO_D
        cg(n)%p%max%hi%k=cg(n-1)%p%max%hi%k/mg_ref_ratio
#else
        cg(n)%p%max%hi%k=1
#endif
        ! use old comms
        cg(n)%t%comm3d = cg(n-1)%t%comm3d
        cg(n)%t%comm = cg(n-1)%t%comm
        cg(n)%t%loc_comm = cg(n-1)%t%loc_comm
        
        cg(n)%t%left = cg(n-1)%t%left
        cg(n)%t%right = cg(n-1)%t%right
        cg(n)%t%top = cg(n-1)%t%top
        cg(n)%t%bottom = cg(n-1)%t%bottom
        cg(n)%t%forward = cg(n-1)%t%forward
        cg(n)%t%behind = cg(n-1)%t%behind
        if (mype==0.and.verbose.gt.4)then
           write(0,*) '[',mype, '] one pe reduction'
           write(0,*) '[',mype, '] X:',cg(n)%t%ipe%i
           write(0,*) '[',mype, '] Y:',cg(n)%t%ipe%j
           write(0,*) '[',mype, '] Z:',cg(n)%t%ipe%k
           write(0,*) '[',mype, '] nx', cg(n)%t%npe%i
           write(0,*) '[',mype, '] ny', cg(n)%t%npe%j
           write(0,*) '[',mype, '] nz', cg(n)%t%npe%k
        endif
     else ! normal split
        cg(n)%t%npe%i = cg(n-1)%t%npe%i/mg_ref_ratio
        cg(n)%t%npe%j = cg(n-1)%t%npe%j/mg_ref_ratio
        cg(n)%t%ipe%i = (cg(n-1)%t%ipe%i-1)/mg_ref_ratio + 1
        cg(n)%t%ipe%j = (cg(n-1)%t%ipe%j-1)/mg_ref_ratio + 1
#ifdef TWO_D
        cg(n)%t%npe%k = cg(n-1)%t%npe%k
        cg(n)%t%ipe%k = cg(n-1)%t%ipe%k
#else
        cg(n)%t%npe%k = cg(n-1)%t%npe%k/mg_ref_ratio
        cg(n)%t%ipe%k = (cg(n-1)%t%ipe%k-1)/mg_ref_ratio + 1
#endif
        ! zero based pe ID, k,j,i
        pe_id = (cg(n)%t%ipe%i-1)*cg(n)%t%npe%j*cg(n)%t%npe%k &
             + (cg(n)%t%ipe%j-1)*cg(n)%t%npe%k + (cg(n)%t%ipe%k-1)
        ii = mod(cg(n-1)%t%ipe%i-1,2)
        jj = mod(cg(n-1)%t%ipe%j-1,2)
        kk = mod(cg(n-1)%t%ipe%k-1,2)
#ifdef TWO_D
        ii = jj + 2*ii           ! local zero based ID
#else 
        ii = kk + 2*jj + 4*ii    ! local zero based ID
#endif
        call MPI_COMM_SPLIT(cg(n-1)%t%comm,ii,pe_id,cg(n)%t%comm,ierr)
        call MPI_COMM_SPLIT(cg(n-1)%t%comm,pe_id,ii,cg(n)%t%loc_comm,ierr)
        
        npes(1) = cg(n)%t%npe%i
        npes(2) = cg(n)%t%npe%j
        npes(3) = cg(n)%t%npe%k
        
        call MPI_cart_create(cg(n)%t%comm,ndims,npes,&
             periodic, .false., cg(n)%t%comm3D, ierr)
        
        ! debug              
        call MPI_comm_rank(cg(n)%t%comm3D,ii,ierr)
        call MPI_Cart_Coords(cg(n)%t%comm3D,ii,ndims,npes,ierr)
        if (cg(n)%t%ipe%i.ne.npes(1)+1) stop '%t%ipe%i'
        if (cg(n)%t%ipe%j.ne.npes(2)+1) stop '%t%ipe%j'
        if (cg(n)%t%ipe%k.ne.npes(3)+1) stop '%t%ipe%k'
        if (mype==mpisize/2.and.verbose.gt.4)then
           write(0,*) '[',mype, '] cart rank',ii
           write(0,*) '[',mype, '] X:',cg(n)%t%ipe%i,npes(1)+1
           write(0,*) '[',mype, '] Y:',cg(n)%t%ipe%j,npes(2)+1
           write(0,*) '[',mype, '] Z:',cg(n)%t%ipe%k,npes(3)+1
           write(0,*) '[',mype, '] npx', cg(n)%t%npe%i
           write(0,*) '[',mype, '] npy', cg(n)%t%npe%j
           write(0,*) '[',mype, '] npz', cg(n)%t%npe%k
        endif

        ! do a box -- todo
        call MPI_Cart_Shift(cg(n)%t%comm3D,0,1,cg(n)%t%left,cg(n)%t%right,ierr)! Neighbors
        call MPI_Cart_Shift(cg(n)%t%comm3D,1,1,cg(n)%t%bottom,cg(n)%t%top,ierr)
        call MPI_cart_Shift(cg(n)%t%comm3D,2,1,cg(n)%t%behind,cg(n)%t%forward,ierr)
        !     logical size of grid stays the same, npe is reduced
        cg(n)%p%max%hi%i=cg(n-1)%p%max%hi%i
        cg(n)%p%max%hi%j=cg(n-1)%p%max%hi%j
        cg(n)%p%max%hi%k=cg(n-1)%p%max%hi%k
        cg(n)%p%max%lo%i=1
        cg(n)%p%max%lo%j=1
        cg(n)%p%max%lo%k=1
     endif
     ! print topology
     if (mype==0.and.verbose.gt.3) then
        write(6,*) '[',mype,'] level ',n,', nxl=',cg(n)%p%max%hi%i,'npx=',cg(n)%t%npe%i
        if (mype==-1) then
           write(6,*) '[',mype,'] ipx=',cg(n)%t%ipe%i, ',ipy=',cg(n)%t%ipe%j,',ipz= ',cg(n)%t%ipe%k
           write(6,*) '[',mype,'] npx=',cg(n)%t%npe%i, ',npy=',cg(n)%t%npe%j,',npz=',cg(n)%t%npe%k
           write(6,*) '[',mype,'] nxl=',&
                cg(n)%p%max%hi%i,',nyl=',&
                cg(n)%p%max%hi%j,  ',nzl=',&
                cg(n)%p%max%hi%k
        endif
     end if
     ncgrids = ncgrids + 1 ! keep track for ease
  enddo ! non-finest coarse grid construction
  !----------------------------------------------------------------------
  ! derived quantities for coarse grids
  do n=ncgrids-1,0,-1
     cg(n)%p%all%lo%i=1-nsg%i
     cg(n)%p%all%lo%j=1-nsg%j
     cg(n)%p%all%hi%i=cg(n)%p%max%hi%i+nsg%i
     cg(n)%p%all%hi%j=cg(n)%p%max%hi%j+nsg%j
#ifndef TWO_D
     cg(n)%p%all%lo%k=1-nsg%k
     cg(n)%p%all%hi%k=cg(n)%p%max%hi%k+nsg%k
#else
     cg(n)%p%all%lo%k=1
     cg(n)%p%all%hi%k=1
#endif 
     if (verbose.gt.1.and.mype==mpisize/2) write (6,'(A,I2,A,I2,I2,I2,A,I4,I4,I4,A,I4,I4,I4,A,I4)'),&
          'new level:',n,&
          ') allocated lo:',cg(n)%p%all%lo%i,cg(n)%p%all%lo%j,cg(n)%p%all%lo%k&
          ,' hi:',cg(n)%p%all%hi%i,cg(n)%p%all%hi%j,cg(n)%p%all%hi%k&
          ,' max:',cg(n)%p%max%hi%i,cg(n)%p%max%hi%j,cg(n)%p%max%hi%k,', npx=',cg(n)%t%npe%i
  end do
  !----------------------------------------------------------------------
  ! setup SR 0 if we have SR
  if (nloc%i*mg_ref_ratio.le.sr_max_loc_sz .and. nloc%j*mg_ref_ratio.le.sr_max_loc_sz &
#ifdef TWO_D
       ) then
#else 
     .and. nloc%k*mg_ref_ratio.le.sr_max_loc_sz) then
#endif
     if (mype==0.and.verbose.gt.3) then
        write(6,*) '[',mype,'] setup SR grid 0, loc n:',nloc%i,nloc%j,nloc%k
     end if
     gsr(0)%t%loc_comm = mpi_comm_null ! 
     gsr(0)%t%comm = comm3d          ! for norms
     gsr(0)%t%comm3d = mpi_comm_null ! flag in BC
     gsr(0)%t%npe%i=NPeAxis(1)       ! for BCs
     gsr(0)%t%npe%j=NPeAxis(2)
     gsr(0)%t%npe%k=NPeAxis(3)
     gsr(0)%t%ipe%i=iPeAxis(1)
     gsr(0)%t%ipe%j=iPeAxis(2)
     gsr(0)%t%ipe%k=iPeAxis(3)
     ! dx
     gsr(0)%p%dx%i = cg(0)%p%dx%i
     gsr(0)%p%dx%j = cg(0)%p%dx%j
#ifdef TWO_D
     gsr(0)%p%dx%k  = 0.d0
#else 
     gsr(0)%p%dx%k = cg(0)%p%dx%k
#endif
     ! new valide size, needs to be grown later
     gsr(0)%p%max%hi%i=nloc%i
     gsr(0)%p%max%hi%j=nloc%j
#ifdef TWO_D
     gsr(0)%p%max%hi%k=1
#else 
     gsr(0)%p%max%hi%k=nloc%k
#endif
     gsr(0)%p%max%lo%i=1
     gsr(0)%p%max%lo%j=1
     gsr(0)%p%max%lo%k=1
  else
     gsr(0)%p%max%hi%i=0 ! flag for empty
     gsr(0)%p%max%hi%j=0
     gsr(0)%p%max%hi%k=0     
  end if
  !----------------------------------------------------------------------
  ! make real SR levels, from coarse to fine
  nsr = 0 ! global var set here
  do n=-1,-dom_max_sr_grids,-1
     ! size of grid, use 'nloc' to keep track of real size
     nloc%i = nloc%i*mg_ref_ratio 
     nloc%j = nloc%j*mg_ref_ratio
     nloc%k = nloc%k*mg_ref_ratio
     gsr(n)%t%loc_comm = mpi_comm_null
     gsr(n)%t%comm = comm3d          ! for norms
     gsr(n)%t%comm3d = mpi_comm_null ! flag in BC
     gsr(n)%t%npe%i=NPeAxis(1)       ! for BCs
     gsr(n)%t%npe%j=NPeAxis(2)
     gsr(n)%t%npe%k=NPeAxis(3)
     gsr(n)%t%ipe%i=iPeAxis(1)
     gsr(n)%t%ipe%j=iPeAxis(2)
     gsr(n)%t%ipe%k=iPeAxis(3)

     ! are we going to make an SR level?
     if (nloc%i.le.sr_max_loc_sz .and. nloc%j.le.sr_max_loc_sz &
#ifdef TWO_D
          ) then
#else 
        .and. nloc%k.le.sr_max_loc_sz) then
#endif
        nsr = nsr + 1 ! have sr
        ! dx
        gsr(n)%p%dx%i = gsr(n+1)%p%dx%i/mg_ref_ratio
        gsr(n)%p%dx%j = gsr(n+1)%p%dx%j/mg_ref_ratio
#ifdef TWO_D
        gsr(n)%p%dx%k  = 0.d0
#else 
        gsr(n)%p%dx%k = gsr(n+1)%p%dx%k/mg_ref_ratio
#endif
        ! new valide size, needs to be grown later
        gsr(n)%p%max%hi%i=nloc%i
        gsr(n)%p%max%hi%j=nloc%j
#ifdef TWO_D
        gsr(n)%p%max%hi%k=1
#else 
        gsr(n)%p%max%hi%k=nloc%k
#endif
        gsr(n)%p%max%lo%i=1
        gsr(n)%p%max%lo%j=1
        gsr(n)%p%max%lo%k=1
     else
        gsr(n)%p%max%hi%i=0 
        gsr(n)%p%max%hi%j=0
        gsr(n)%p%max%hi%k=0
     end if
  end do
  if (mype==0.and.verbose.gt.3.and.nsr.gt.0) then
     write(6,*) '[',mype,'] setup ',nsr,' SR levels, buffer schedual:',sr_base_bufsz,sr_bufsz_inc
  end if
  !----------------------------------------------------------------------
  ! add buffer/ghosts, set size, from fine to coarse
  if (nsr.gt.0) then
     srbf = sr_base_bufsz-sr_bufsz_inc ! the buffer schedual
     do n=-nsr,0,1
        srbf = srbf+sr_bufsz_inc ! number of buffer cells
        if (srbf.gt.gsr(n)%p%max%hi%i.or.srbf.gt.gsr(n)%p%max%hi%j&
#ifndef TWO_D
             .or.srbf.gt.gsr(n)%p%max%hi%k)&
#endif
             stop 'SR buffer larger than grid'
        ! offset to line up with finer grid for R & P: c.bfsz - f.bfsz/2 = c.bfsz/2 + inc/2
        gsr(n)%cfoffset%i = srbf - (srbf-sr_bufsz_inc)/2 
        gsr(n)%cfoffset%j = srbf - (srbf-sr_bufsz_inc)/2
        gsr(n)%cfoffset%k = srbf - (srbf-sr_bufsz_inc)/2
        
        ! clip buffers at BCs
        ii = srbf
        jj = srbf
#ifndef XPERIODIC
        if (iPeAxis(1)==1) then
           ii = 0
           gsr(n)%cfoffset%i = 0
        endif
        if (iPeAxis(1)==NPeAxis(1)) then
           jj = 0
        endif
#endif
        ! val: valid region inside normal grid
        gsr(n)%val%lo%i=ii+1
        gsr(n)%val%hi%i=gsr(n)%p%max%hi%i+ii      ! size of valid + low buffer size
        gsr(n)%p%max%hi%i=gsr(n)%p%max%hi%i+ii+jj ! new comp area
        ! Y
        ii = srbf
        jj = srbf
#ifndef YPERIODIC
        if (iPeAxis(2)==1) then
           ii = 0
           gsr(n)%cfoffset%j = 0
        endif
        if (iPeAxis(2)==NPeAxis(2)) then
           jj = 0
        endif
#endif
        gsr(n)%val%lo%j=ii+1
        gsr(n)%val%hi%j=gsr(n)%p%max%hi%j+ii
        gsr(n)%p%max%hi%j=gsr(n)%p%max%hi%j+ii+jj ! new comp area
        ! Z
#ifdef TWO_D
        gsr(n)%p%max%hi%k=1
        gsr(n)%val%lo%k=1
        gsr(n)%val%hi%k=1
#else
        ii = srbf
        jj = srbf
#ifndef ZPERIODIC
        if (iPeAxis(3)==1) then
           ii = 0
           gsr(n)%cfoffset%k = 0
        endif
        if (iPeAxis(3)==NPeAxis(3)) then
           jj = 0
        endif
#endif
        gsr(n)%val%lo%k=ii+1
        gsr(n)%val%hi%k=gsr(n)%p%max%hi%k+ii
        gsr(n)%p%max%hi%k=gsr(n)%p%max%hi%k+ii+jj
#endif 
        ! data size
        gsr(n)%p%all%lo%i=-nsg%i+1
        gsr(n)%p%all%hi%i= gsr(n)%p%max%hi%i+nsg%i
        gsr(n)%p%all%lo%j=-nsg%j+1
        gsr(n)%p%all%hi%j= gsr(n)%p%max%hi%j+nsg%j
#ifdef TWO_D
        gsr(n)%p%all%lo%k=1
        gsr(n)%p%all%hi%k=1
#else
        gsr(n)%p%all%lo%k=-nsg%k+1
        gsr(n)%p%all%hi%k=gsr(n)%p%max%hi%k+nsg%k
#endif
        if (n==-nsr) then ! no cf needed for finest
           gsr(n)%cfoffset%i = 0
           gsr(n)%cfoffset%j = 0
           gsr(n)%cfoffset%k = 0
        end if
     end do
  end if
  ! print SR levels
  if (nsr>0.and.verbose.gt.1.and.mype==mpisize/2) then
     do n=0,-nsr,-1
        write (6,'(A,I2,A,I3,I3,I3,A,I3,I3,I3,A,I4,I4,I4,A,I3,I3,I3)'),&
             'new SR level:',n,') valid lo:',gsr(n)%val%lo%i,gsr(n)%val%lo%j,gsr(n)%val%lo%k&
             ,' hi:',gsr(n)%val%hi%i,gsr(n)%val%hi%j,gsr(n)%val%hi%k&
             ,' max:',gsr(n)%p%max%hi%i,gsr(n)%p%max%hi%j,gsr(n)%p%max%hi%k,&
             ', (crs) R/P offset:',gsr(n)%cfoffset%i,gsr(n)%cfoffset%j,gsr(n)%cfoffset%k
     end do
  end if
  return
end subroutine new_grids_private
