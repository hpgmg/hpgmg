!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
! domain size & BCs
module domain
  save
  double precision:: xl,xr,yl,yr,zl,zr
  logical:: periodic(3)
end module domain
!  topology and solver and general parametersdestroy
module pms
  ! 
  integer :: problemType
  ! topo
  integer,parameter:: dom_max_grids=20
  integer,parameter:: dom_max_sr_grids=10
  integer,parameter:: mg_ref_ratio=2
  integer:: bot_min_sz ! min size for bottom solver
  integer:: pe_min_sz  ! min size box to start aggregating
  integer:: sr_max_loc_sz    ! min size of box for finest grid in SR coarse grid solver
  integer:: nsr ! number of SR grids, used to size grid array (0 is finest 'sr coarse grid')
  integer:: ncgrids ! number of grids in SR coarse grid (0:ncgrids-1)
  integer:: sr_base_bufsz,sr_bufsz_inc ! the SR buffer schedual
  ! solver
  integer:: nvcycles ! 0 for pure FMG
  integer:: nfcycles ! 0 for no FMG
  integer:: ncycles  ! V==1, W==2
  integer:: nfmgvcycles ! 1 for normal FMG, more for more power
  double precision:: rtol
  integer:: nsmoothsup 
  integer:: nsmoothsdown
  integer:: ncoarsesolveits
  ! general
  integer:: verbose 
  integer:: num_solves
end module pms
!-----------------------------------------------------------------------
module iounits
  integer, parameter:: irun=1
  integer, parameter:: ibinoutput=2
  integer, parameter:: iconv=3
end module iounits
!-----------------------------------------------------------------------
module tags
  integer:: MSG_XCH_XLOW_TAG,MSG_XCH_XHI_TAG
  integer:: MSG_XCH_YLOW_TAG,MSG_XCH_YHI_TAG
  integer:: MSG_XCH_ZLOW_TAG,MSG_XCH_ZHI_TAG
end module tags
!-----------------------------------------------------------------------
module mpistuff
#ifdef HAVE_PETSC
#include "finclude/petscdef.h" 
#include "finclude/petsclogdef.h"  
  use petsc
  double precision :: flops
  integer :: events(10)
#else
  include "mpif.h"
#endif
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
module pe_patch_data_module
  type pe_patch
     integer:: imax,jmax,kmax                       ! compute region w/0 ghosts
     integer:: ivallo,jvallo,kvallo,ivalhi,jvalhi,kvalhi ! valid region
     integer:: ilo,ihi,jlo,jhi,klo,khi              ! data size w/ ghosts
     integer:: left,right,top,bottom,behind,forward ! cache of pe neighbors
     integer:: comm3d,loc_comm,comm ! 3d comm, local comm of redunent parents, comm (split)
     integer:: ipex,ipey,ipez ! my pe in 3D index space
     integer:: npex,npey,npez ! size if pe index space
     double precision:: dxg,dyg,dzg ! why 'g'?
  end type pe_patch
end module pe_patch_data_module
!-----------------------------------------------------------------------
module GridModule
  use pe_patch_data_module
  implicit none
  ! PDE, Disc
  integer,parameter:: nsg=1  ! number of stencil ghosts     
  integer,parameter:: nvar=1 ! num vars per cell
contains
  !-----------------------------------------------------------------
  logical function is_top(g)
    use pms, only:bot_min_sz
    implicit none
    type(pe_patch):: g
    is_top = ( &
         (g%imax.le.bot_min_sz) .or. &
         (g%jmax.le.bot_min_sz) &
#ifndef TWO_D
         .or. (g%kmax.le.bot_min_sz) &
#endif
         )
  end function is_top
  !-----------------------------------------------------------------
  integer function getIglobalx(g)
    use mpistuff, only:mype
    implicit none
    type(pe_patch):: g
    getiglobalx = (g%ipex-1)*g%imax + 1 ! one based of my first index
  end function getIglobalx
 !-----------------------------------------------------------------
  integer function getIglobaly(g)
    implicit none
    type(pe_patch):: g
    getiglobaly = (g%ipey-1)*g%jmax + 1 
  end function getIglobaly
 !-----------------------------------------------------------------
  integer function getIglobalz(g)
    implicit none
    type(pe_patch):: g
    getiglobalz = (g%ipez-1)*g%kmax + 1 
  end function getIglobalz
  !-----------------------------------------------------------------
  subroutine new_grids(g,NPeAxis,iPeAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)
    use pms, only:dom_max_sr_grids,dom_max_grids
    implicit none    
    integer,intent(in):: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
    integer :: NPeAxis(3),iPeAxis(3)
    type(pe_patch),intent(out):: g(-dom_max_sr_grids:dom_max_grids-1)
    interface       
       subroutine new_grids_private(g,NPeAxis,iPeAxis,nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d)
         use pms, only:dom_max_sr_grids,dom_max_grids
         use pe_patch_data_module
         implicit none
         integer,intent(in):: nx,ny,nz,nxlocal,nylocal,nzlocal,comm3d
         integer :: NPeAxis(3),iPeAxis(3)
         type(pe_patch),intent(out):: g(-dom_max_sr_grids:dom_max_grids-1) ! this needs to be set but max_grids is not visable!!  
       end subroutine new_grids_private
    end interface

    call new_grids_private(g,NPeAxis,iPeAxis,nx,ny,nz,&
         nxlocal,nylocal,nzlocal,comm3d)

  end subroutine new_grids
  !-----------------------------------------------------------------
  subroutine destroy_grids(g)
    use pms, only:dom_max_sr_grids,dom_max_grids
    use pe_patch_data_module
    implicit none    
    type(pe_patch):: g(-dom_max_sr_grids:dom_max_grids-1)
    interface       
       subroutine destroy_grids_private(g)
         use pms, only:dom_max_sr_grids,dom_max_grids
         use pe_patch_data_module
         implicit none
         type(pe_patch):: g(-dom_max_sr_grids:dom_max_grids-1)
       end subroutine destroy_grids_private
    end interface
    
    call destroy_grids_private(g)
    
  end subroutine destroy_grids
  !
end module GridModule

