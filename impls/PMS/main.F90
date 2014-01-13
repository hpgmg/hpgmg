!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
program main
  use GridModule
  use domain, only:periodic
  use pms, only:verbose,dom_max_grids,dom_max_sr_grids
  use mpistuff
  implicit none
  integer:: rank,comm3d,ndims
  integer:: NPeAxis(3),iPeAxis(3),npex,npey,npez,nx,ny,nz,nxlocal,nylocal,nzlocal
  type(pe_patch)::grids(dom_max_grids+dom_max_sr_grids)

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

  ! get fine grid parameters, problem type, and global variables
  call get_params(npe, npex, npey, npez, nx, ny, nz, &
       nxlocal, nylocal, nzlocal)

  ! Initialize the fine grid communicator, set periodicity
  ndims = 3
  NPeAxis(1) = npex
  NPeAxis(2) = npey
  NPeAxis(3) = npez
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
  if((mype==0.and.verbose.gt.3).or.mype.ne.rank) then
     if (verbose.gt.0) write(6,*) 'main: warning old mype=',mype,', comm3D new rank = ',rank
  endif
  mype=rank

  call MPI_Cart_Coords(comm3d,rank,ndims,iPeAxis,ierr)
  iPeAxis(1) = iPeAxis(1) + 1 ! one based
  iPeAxis(2) = iPeAxis(2) + 1
  iPeAxis(3) = iPeAxis(3) + 1
  
  ! run it
  call SetupDomain() ! sets size of domain
  
  ! sets global vars
  call new_grids(grids,NPeAxis,iPeAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)
  
  ! do it 
  call driver(grids)
  
  ! end it
  call destroy_grids_private(grids)
  call FinalizeIO()
  call MPI_Finalize(ierr)
end program MAIN
!-----------------------------------------------------------------
!  driver - the PETSc main
!-----------------------------------------------------------------
subroutine driver(grids)
  use iounits
  use mpistuff 
  use pms
  use GridModule
  use error_data_module
  implicit none
  !
  type(pe_patch)::grids(-dom_max_sr_grids:dom_max_grids-1)
  !
  integer:: isolve,ii,coarsest_grid,nViters,kk
  integer,parameter:: dummy_size=1024*1024 ! cache flush array size
  double precision:: dummy(dummy_size),res0,rateu,rategrad
  type(error_data)::errors(-nsr:ncgrids-1+nvcycles)
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
                 if (verbose .gt.0) write(6,'(A,I2,A,E14.6,A,E14.6)') &
                      'driver: level ',ii, ': converg order u=', &
                      log(rateu)*log2r,', grad(u)=',log(rategrad)*log2r
                 write(iconv,900) isolve,log(rateu)*log2r,log(rategrad)*log2r,errors(ii)%uerror,&
                      errors(ii)%graduerr,errors(ii)%resid
                 
                 if( is_top(grids(ii)) ) then 
                    coarsest_grid = ii
                    exit ! grids and errors are in reverse order - this is just a counter
                 end if
              end do
           end if
 
           ! print V cycles
           do ii=coarsest_grid,coarsest_grid+nViters-1 ! i is zero bases for 'errors'
              if (verbose.gt.0) write(6,'(A,I2,A,E14.6,A,E14.6)') 'driver: V cycle',&
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
#ifdef HAVE_PETSC
  call PetscFinalize(ierr)
#endif
  return
end subroutine driver
!-----------------------------------------------------------------
! get command line args, set problemType and fine grid sizes
!-----------------------------------------------------------------       
subroutine get_params(npes, npex, npey, npez, &
     nx, ny, nz, nxlocal, nylocal, nzlocal)
  use pms
  use mpistuff
  use tags
  !use tags, only:msg_tag_inc
#ifndef HAVE_COMM_ARGS
  use f2kcli        ! command line args class
#endif
  implicit none
  integer,intent(out) :: npex,npey,npez,nx,ny,nz,nxlocal,nylocal,nzlocal
  integer,intent(in) :: npes
  !       
  CHARACTER(LEN=256) :: LINE
  CHARACTER(LEN=32)  :: CMD,CMD2
  INTEGER            :: NARG,IARG,fact,nsmooths
  integer*8 :: ti
  integer, parameter :: nargs=22
  double precision t2,t1,iarray(nargs)
  ! default input args: ?local, n?, ?pes
  npex=-1
  npey=-1
  npez=-1
  nx=-1
  ny=-1
  nz=-1
  nxlocal=-1
  nylocal=-1
  nzlocal=-1
  !       Problem type 
  !       0: (x^4-x^2)^D - Diri BCs
  !       1: bubble      - periodic domain???
  problemType = 0
  ! solver parameters
  nvcycles = 0 ! pure FMG
  nfcycles = 1 ! pure FMG
  rtol = 1.d-5    ! only used for V-cycles
  ncycles = 1     ! v-cycles or w-cycles
  nfmgvcycles = 1 ! no interface for this (always 1)
  nsmooths = 2
  verbose = 1
  bot_min_sz = 2 ! min size for bottom solver (grad get messed up with 1)
  pe_min_sz = 8  ! or 16 ...
  num_solves = 1
  sr_max_loc_sz = 0 ! no SR by default
  sr_base_bufsz = 3
  sr_bufsz_inc = 1
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
        
        if( CMD == '-nxpe' ) READ (CMD2, '(I10)') npex
        if( CMD == '-nype' ) READ (CMD2, '(I10)') npey
        if( CMD == '-nzpe' ) READ (CMD2, '(I10)') npez

        if( CMD == '-bot_min_sz' ) READ (CMD2, '(I10)') bot_min_sz
        if( CMD == 'pe_min_sz' ) READ (CMD2, '(I10)') pe_min_sz

        if( CMD == '-problem' ) READ (CMD2, '(I10)') problemType

        if( CMD == '-nvcycles' ) READ (CMD2, '(I10)') nvcycles
        if( CMD == '-nsmooths' ) READ (CMD2, '(I10)') nsmooths
        if( CMD == '-nfcycles' ) READ (CMD2, '(I10)') nfcycles
        if( CMD == '-rtol' )    READ (CMD2, '(E14.6)') rtol
        if( CMD == '-ncycle' )  READ (CMD2, '(I10)') ncycles

        if( CMD == '-verbose' )  READ (CMD2, '(I10)') verbose
        if( CMD == '-num_solves' )  READ (CMD2, '(I10)') num_solves

        if( CMD == '-sr_max_loc_sz' )  READ (CMD2, '(I10)') sr_max_loc_sz

        if( CMD == '-sr_base_bufsz' )  READ (CMD2, '(I10)') sr_base_bufsz
        if( CMD == '-sr_bufsz_inc' )  READ (CMD2, '(I10)') sr_bufsz_inc
     END DO
#ifdef TWO_D
     nz = 1; nzlocal = 1; npez = 1;
#endif
     iarray(1) = nx
     iarray(2) = ny
     iarray(3) = nz
     
     iarray(4) = nxlocal
     iarray(5) = nylocal
     iarray(6) = nzlocal
     
     iarray(7) = npex
     iarray(8) = npey
     iarray(9) = npez
     
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
     
     npex = int(iarray(7))
     npey = int(iarray(8))
     npez = int(iarray(9))
     
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
  endif
  if (nfcycles .ne. 0) nfcycles = 1 ! only valid 0,1
  if (nvcycles .lt. 0) nvcycles = 0 ! only valid 0,1
  if (nfcycles+nvcycles .le. 0) print *, 'no solver iterations!!!'
  ! general constructor stuff
  nsmoothsup = nsmooths
  nsmoothsdown = nsmooths
  ncoarsesolveits = bot_min_sz*bot_min_sz
  MSG_XCH_XLOW_TAG=100
  MSG_XCH_XHI_TAG=200
  MSG_XCH_YLOW_TAG=300
  MSG_XCH_YHI_TAG=400
  MSG_XCH_ZLOW_TAG=500
  MSG_XCH_ZHI_TAG=600

  ! get ijk pes
  t1 = npes
#ifdef TWO_D
  npez = 1
  nz = 1
  nzlocal = 1
  if (npex == -1) then
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
     npex=fact*floor(dsqrt(t1/fact))
     npey=t1/npex
     if (npes .ne. npez*npex*npey) then
        write(6,*) 'INCORRECT NP:', NPES,'npex=',npex,&
             'npey=',npey,'npez=',npez,'npes=',npes
        stop
     endif
  else (npey == -1) then
     npey=t1/npex
     if(npey<1)stop'too many npex???'
     fact = npex/npey
     if( fact < 1 ) then
        write(6,*) 'npex must be greater than npey', NPES,&
       'npex=',npex,&
       'npey=',npey,'npez=',npez,'npes=',npes
        stop
     endif
  endif
#else
  if (npex == -1) then
     t1 = npes
     if(npes==1) then ! one pe
        npex = 1; npey = 1; npez = 1;
     elseif(nx==ny .and. ny==nz)then ! square
        npex = floor(t1**.3333333334d0)
        npey = npex
        npez = npes/(npex*npey)
     elseif(nz==ny)then ! z==y
        t1 = ny*npes; t2 = nx; t1 = t1/t2 
        npey = floor(t1**.333333334d0)
        npez = npey
        npex = npes/(npey*npez) 
     elseif(nz==nx)then
        t1 = nz*npes; t2 = ny; t1 = t1/t2 
        npex = floor(t1**.333333334d0)
        npez = npex
        npey = npes/(npex*npez) 
     elseif(nx==ny)then
        t1 = ny*npes; t2 = nz; t1 = t1/t2 
        npey = floor(t1**.333333334d0)
        npex = npey
        npez = npes/(npey*npex) 
     else
        stop 'need to specify x/y/z pe'
     endif
  endif
#endif     
  if (npes .ne. npez*npex*npey) then
     write(6,*) 'INCORRECT NP: npex=',npex, &
          'npey=',npey,'npez=',npez,'npes=',npes
     stop
  endif
  ! x
  if (nxlocal.eq.-1 .and. nx.eq.-1) nxlocal = 64        ! default
  if (nxlocal.eq.-1 .and. nx.ne.-1) nxlocal = nx/npex 
  if (nxlocal.ne.-1 .and. nx.eq.-1) nx = nxlocal*npex 
  if (nx .ne. nxlocal*npex) stop 'nx .ne. nxlocal*npex'
  ! y
  if (nylocal.eq.-1 .and. ny.eq.-1) nylocal = nxlocal   ! default square
  if (nylocal.eq.-1 .and. ny.ne.-1) nylocal = ny/npey 
  if (nylocal.ne.-1 .and. ny.eq.-1) ny = nylocal*npey 
  if (ny .ne. nylocal*npey) stop 'ny .ne. nylocal*npey'
  ! z
  if (nzlocal.eq.-1 .and. nz.eq.-1) nzlocal = nxlocal   ! default square
  if (nzlocal.eq.-1 .and. nz.ne.-1) nzlocal = nz/npez 
  if (nzlocal.ne.-1 .and. nz.eq.-1) nz = nzlocal*npez 
  if (nz .ne. nzlocal*npez) stop 'nz .ne. nzlocal*npez'
  ! print topology
  if(mype==0 .and. verbose .gt. 3) then
     write(6,*) '[',mype,'] nx =',nx,     ',ny= ',ny     ,',nz= ',nz
     write(6,*) '[',mype,'] npx=',npex, ',npy=',npey, ',npz=',npez
     write(6,*) '[',mype,'] nxl=',nxlocal,',nyl=',nylocal,',nzl=',nzlocal
  endif
end subroutine get_params

!-----------------------------------------------------------------
subroutine SetupDomain()
  use domain
  use pms, only:problemType
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
