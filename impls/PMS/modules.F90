!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
module domain
  save
  double precision:: xl,xr,yl,yr,zl,zr
  integer :: problemType
  !  double precision:: dx,dy,dz
  !  double precision, allocatable,dimension(:):: xc,yc,zc
end module domain
!-----------------------------------------------------------------------
module iounits
  integer, parameter:: irun=1
  integer, parameter:: ibinoutput=2
  integer, parameter:: iconv=3
end module iounits
!-----------------------------------------------------------------------
module tags
  integer,parameter:: MSG_XCH_XLOW_TAG=1,MSG_XCH_XHI_TAG=2
  integer,parameter:: MSG_XCH_YLOW_TAG=3,MSG_XCH_YHI_TAG=4
  integer,parameter:: MSG_XCH_ZLOW_TAG=5,MSG_XCH_ZHI_TAG=6
  integer,parameter:: MSG_MAX_TAG=10
end module tags
!-----------------------------------------------------------------------
module mpistuff
  include "mpif.h"
  integer:: status(MPI_STATUS_SIZE),ierr
  integer:: mype ,npe
  integer,parameter::ERROR_CARTCOORDS=1
  integer,parameter::ERROR_CARTSHIFT=2
  integer,parameter::ERROR_WAIT=3
  integer,parameter::ERROR_SEND=4
  integer,parameter::ERROR_RECV=5
  integer,parameter::ERROR_ALLREDUCE=6
end module mpistuff
!-----------------------------------------------------------------------
! error data module
module error_data_module
  type error_data
     double precision:: uerror,graduerr,resid
  end type error_data
  integer :: err_lev ! cache of currant level for error output
  integer,parameter::error_norm=3 ! 3 == inf
end module error_data_module
!-----------------------------------------------------------------------
! Grid module
module proc_patch_data_module
  type proc_patch
     integer:: imax,jmax,kmax                       ! size w/0 ghosts
     integer:: ilo,ihi,jlo,jhi,klo,khi              ! data size w/ ghosts
     integer:: left,right,top,bottom,behind,forward ! cache of proc neighbors
     integer:: comm3d,comm,loc_comm ! 3d comm, base (mpi_comm_world or split), comm with r^D pes
     integer:: iprocx,iprocy,iprocz ! my proc in 3D index space
     integer:: nprocx,nprocy,nprocz ! size if proc index space
     integer:: iglobalx,iglobaly,iglobalz ! easy to generate, could remove
     double precision:: dxg,dyg,dzg ! why 'g'?
  end type proc_patch
end module proc_patch_data_module
!-----------------------------------------------------------------------
#define MAX_GRIDS 20
module GridModule
  use proc_patch_data_module
  implicit none
  ! PDE,Disc
  integer,parameter:: ng=1     
  integer,parameter:: nvar=1   
  ! topo
  integer,parameter:: max_grids=MAX_GRIDS
  integer:: bot_min_size ! min size for bottom solver
  integer:: mg_min_size  ! min size box to start aggregating
  integer:: ncoarsesolveits
  ! solver
  integer:: nvcycles ! 0 for pure FMG
  integer:: nfcycles ! 0 for no FMG
  integer:: ncycles  ! V==1, W==2
  integer:: nfmgvcycles ! 1 for normal FMG, more for more power
  double precision:: rtol
  integer:: nsmoothsup 
  integer:: nsmoothsdown
  ! general
  integer:: verbose 
contains
  !-----------------------------------------------------------------
  logical function is_top(g)
    implicit none
    type(proc_patch):: g
    is_top = ( &
         (g%imax.le.bot_min_size .and. g%nprocx.eq.1) .or. &
         (g%jmax.le.bot_min_size .and. g%nprocy.eq.1) &
#ifndef TWO_D
         .or. (g%kmax.le.bot_min_size .and. g%nprocz.eq.1) &
#endif
         )
  end function is_top
  !-----------------------------------------------------------------
  integer function getIglobalx(g)
    use mpistuff, only:mype
    implicit none
    type(proc_patch):: g
    getiglobalx = (g%iprocx-1)*g%imax + 1 ! one based of my first index
  end function getIglobalx
 !-----------------------------------------------------------------
  integer function getIglobaly(g)
    implicit none
    type(proc_patch):: g
    getiglobaly = (g%iprocy-1)*g%jmax + 1 
  end function getIglobaly
 !-----------------------------------------------------------------
  integer function getIglobalz(g)
    implicit none
    type(proc_patch):: g
    getiglobalz = (g%iprocz-1)*g%kmax + 1 
  end function getIglobalz
  !-----------------------------------------------------------------
  subroutine new_grids(g,NProcAxis,iProcAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)
    implicit none    
    integer,intent(in):: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
    integer :: NProcAxis(3),iProcAxis(3)
    type(proc_patch),intent(out):: g(0:max_grids-1)
    interface       
       subroutine new_grids_private(g,NProcAxis,iProcAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)
         use proc_patch_data_module
         implicit none
         integer,intent(in):: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
         integer :: NProcAxis(3),iProcAxis(3)
         type(proc_patch),intent(out):: g(0:MAX_GRIDS-1) ! this needs to be set but max_grids is not visable!!  
       end subroutine new_grids_private
    end interface

    call new_grids_private(g,NProcAxis,iProcAxis,nx,ny,nz,&
         nxlocal,nylocal,nzlocal,comm3d)

  end subroutine new_grids
  !-----------------------------------------------------------------
  subroutine destroy_grids(g,max_sz)
    use proc_patch_data_module
    implicit none    
    integer,intent(in):: max_sz
    type(proc_patch):: g(0:max_sz-1)
    interface       
       subroutine destroy_grids_private(g,max_sz)
         use proc_patch_data_module
         implicit none
         integer,intent(in):: max_sz
         type(proc_patch):: g(0:max_sz-1)       
       end subroutine destroy_grids_private
    end interface
    
    call destroy_grids_private(g,max_sz)
    
  end subroutine destroy_grids
  !
end module GridModule

