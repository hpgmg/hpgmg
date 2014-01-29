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
!  topology, solver, and general parametersdestroy
module pms
  ! 
  integer :: problemType
  ! topo
  integer,parameter:: dom_max_grids=20
  integer,parameter:: dom_max_sr_grids=10
  integer,parameter:: mg_ref_ratio=2
  integer:: bot_min_sz ! min size for bottom solver
  integer:: pe_min_sz  ! min size box to start aggregating
  integer:: sr_max_loc_sz ! max size of box for finest grid in SR
  integer:: nsr ! number of SR grids, used to size grid array (0 is finest 'sr coarse grid')
  integer:: ncgrids ! number of grids in SR coarse grid (0:ncgrids-1)
  integer:: sr_base_bufsz,sr_bufsz_inc ! the SR buffer schedual
  ! solver
  integer:: nvcycles ! 0 for pure FMG
  integer:: nfcycles ! 0 for no FMG
  integer:: ncycles  ! V==1, W==2
  integer:: nfmgvcycles ! 1 for normal FMG, more for more power
  double precision:: rtol
  integer::nsmoothsup 
  integer::nsmoothsdown
  integer::nsmoothsfmg
  integer::ncoarsesolveits
  ! general
  integer:: verbose 
  integer:: num_solves
end module pms
!-----------------------------------------------------------------------
module iounits
  integer, parameter:: istderr=0
  integer, parameter:: istdin=5
  integer, parameter:: istdout=6
  integer, parameter:: irun=8
  integer, parameter:: imesh=16
  integer, parameter:: ifluid=21
  integer, parameter:: idump=71
  integer, parameter:: ibinoutput=41
  integer, parameter:: itecoutput=61
  integer, parameter:: itime=81
  integer, parameter:: idiag=91
  integer, parameter:: iconv=92
  integer, parameter:: imass=93
  integer, parameter:: iEng=94
  integer, parameter:: ipflux=95
  integer, parameter:: iflux=97
  integer, parameter:: iprof=999
  integer, parameter:: itflux=96
end module iounits
!-----------------------------------------------------------------------
module tags
  integer:: MSG_XCH_XLOW_TAG,MSG_XCH_XHI_TAG
  integer:: MSG_XCH_YLOW_TAG,MSG_XCH_YHI_TAG
  integer:: MSG_XCH_ZLOW_TAG,MSG_XCH_ZHI_TAG
  integer:: MSG_RESTRICT_TAG
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
  integer::mype,mpisize
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
     double precision:: uerror,resid
  end type error_data
  integer::err_lev ! cache of currant level for error output
  integer::error_norm ! 3 == inf
end module error_data_module
!-----------------------------------------------------------------------
! base objects
module base_data_module
  !
  type ipoint
     integer:: i,j,k
  end type ipoint
  !
  type dpoint
     double precision:: i,j,k
  end type dpoint
  !  
  type box
     type(ipoint):: lo,hi
  end type box
  ! mpi meta data for level
  type topot
     type(ipoint)::ipe ! my pe in 3D index space
     type(ipoint)::npe ! size if pe index space
     integer::left,right,top,bottom,behind,forward ! my neighbors - should be int[-1:1]^3 
     integer::comm3d   ! 3d comm
     integer::loc_comm ! local comm of redunent parents
     integer::comm     ! comm (split) & norms
  end type topot
  ! meta data for operators (apply,relax)
  type patcht
     type(box)::max     ! compute region w/0 ghosts, lo=(1,1,1) now
     type(box)::all     ! data size w/ ghosts
     type(dpoint)::dx 
  end type patcht
  ! normal patch: patch meta data and processor topo stuff, could split this up more
  type crs_patcht
     type(patcht)::p  ! local size,dx,compute region
     type(topot)::t   ! MPI stuff
  end type crs_patcht
  !
  type sr_patcht
     type(patcht)::p   ! local size,dx,compute region
     type(topot)::t    ! no need for MPI topo in (pure) SR but need ipe for BCs and a comm for norms
     type(ipoint)::cfoffset ! offset to line up with finer grid for R & P
     type(box)::val    ! valid region, norms, zero level copy ...
  end type sr_patcht
  ! array pointer for making arrays of array pointers
  type data_ptr
     double precision,DIMENSION(:,:,:,:),POINTER::ptr
  end type data_ptr
end module base_data_module
!-----------------------------------------------------------------------
module discretization
  use base_data_module
  ! number of stencil ghosts, need to hack for one SR exch.     
  type(ipoint),parameter::nsg=ipoint(1,1,1) 
  ! num vars per cell
  integer,parameter::nvar=1
  ! physics
  double precision,parameter::alpha=-1.d0
end module discretization
!-----------------------------------------------------------------------
! patch data types
module bc_module
  use base_data_module
  interface
     subroutine SetBCs(ux,p,t,ng)
       use discretization,only:nvar
       use base_data_module
       implicit none
       type(ipoint),intent(in)::ng
       type(topot),intent(in)::t
       type(patcht),intent(in)::p
       double precision:: ux(&
            p%all%lo%i:p%all%hi%i,&
            p%all%lo%j:p%all%hi%j,&
            p%all%lo%k:p%all%hi%k,nvar)
     end subroutine SetBCs
     subroutine sr_exchange(u_sr,srg0,cg0)
       use mpistuff,only:mype
       use discretization
       implicit none
       type(sr_patcht),intent(in)::srg0
       type(crs_patcht),intent(in)::cg0
       double precision,dimension(&
            srg0%p%all%lo%i:srg0%p%all%hi%i,&
            srg0%p%all%lo%j:srg0%p%all%hi%j,&
            srg0%p%all%lo%k:srg0%p%all%hi%k,nvar)::u_sr
     end subroutine sr_exchange
  end interface
  type(ipoint)::bc_valid
contains
  !-----------------------------------------------------------------
  integer function getIglobalx(val,ip)
    implicit none
    type(box),intent(in)::val
    type(ipoint),intent(in)::ip
    getiglobalx = (ip%i-1)*(val%hi%i-val%lo%i+1) + 1 ! one based of my first index
  end function getIglobalx
 !-----------------------------------------------------------------
  integer function getIglobaly(val,ip)
    implicit none
    type(box),intent(in)::val
    type(ipoint),intent(in)::ip
    getiglobaly = (ip%j-1)*(val%hi%j-val%lo%j+1) + 1 
  end function getIglobaly
 !-----------------------------------------------------------------
  integer function getIglobalz(val,ip)
    implicit none
    type(box),intent(in)::val
    type(ipoint),intent(in)::ip
    getiglobalz = (ip%k-1)*(val%hi%k-val%lo%k+1) + 1 
  end function getIglobalz
end module bc_module
!-----------------------------------------------------------------------
module grid_module
  use bc_module
  ! interfaces - kernels
  interface
     ! prolongate
     subroutine Prolong(uxF,uxC,ft,ct,fp,cp,cOffset,high)
       !  uxF = uxF + P * uxC
       use discretization
       implicit none
       type(patcht),intent(in)::fp,cp
       type(ipoint),intent(in)::cOffset
       type(topot),intent(in)::ct,ft
       logical,intent(in)::high
       double precision,intent(in)::uxC(&
            cp%all%lo%i:cp%all%hi%i,&
            cp%all%lo%j:cp%all%hi%j,&
            cp%all%lo%k:cp%all%hi%k,nvar)
       double precision,intent(out)::uxF(&
            fp%all%lo%i:fp%all%hi%i,&
            fp%all%lo%j:fp%all%hi%j,&
            fp%all%lo%k:fp%all%hi%k,nvar)
     end subroutine Prolong
     ! restriction
     subroutine RestrictFuse(cp,fp,ct,cOffset,uC,uC2,ux,ux2)
       use discretization
       implicit none
       type(patcht),intent(in)::fp,cp
       double precision,intent(in),dimension(&
            fp%all%lo%i:fp%all%hi%i,&
            fp%all%lo%j:fp%all%hi%j,&
       fp%all%lo%k:fp%all%hi%k,nvar)::ux,ux2
       double precision,intent(out),dimension(&
       cp%all%lo%i:cp%all%hi%i,&
       cp%all%lo%j:cp%all%hi%j,&
       cp%all%lo%k:cp%all%hi%k,nvar)::uC,uC2
       type(ipoint),intent(in)::cOffset
       type(topot),intent(in)::ct
     end subroutine RestrictFuse
     !
     subroutine Apply_const_Lap(uxo,ux,p,t)
       use discretization
       implicit none
       type(patcht)::p
       double precision,intent(out):: uxo(&
            p%all%lo%i:p%all%hi%i,&
            p%all%lo%j:p%all%hi%j,&
            p%all%lo%k:p%all%hi%k,nvar)
       double precision,intent(in)::ux(&
            p%all%lo%i:p%all%hi%i,&
            p%all%lo%j:p%all%hi%j,&
            p%all%lo%k:p%all%hi%k,nvar)
       type(topot),intent(in)::t
     end subroutine Apply_const_Lap
     !
     subroutine GS_RB_const_Lap(phi,rhs,p,t,nits)
       use discretization
       use pms, only:verbose
       implicit none
       type(patcht)::p
       double precision,dimension(&
            p%all%lo%i:p%all%hi%i,&
            p%all%lo%j:p%all%hi%j,&
            p%all%lo%k:p%all%hi%k,nvar)::phi,rhs
       type(topot),intent(in)::t
       integer,intent(in):: nits
     end subroutine GS_RB_const_Lap
     !
     subroutine GS_Lex_const_Lap(phi,rhs,p,t,nits)
       use discretization
       use pms, only:verbose
       implicit none
       type(patcht)::p
       double precision,dimension(&
            p%all%lo%i:p%all%hi%i,&
            p%all%lo%j:p%all%hi%j,&
            p%all%lo%k:p%all%hi%k,nvar)::phi,rhs
       type(topot),intent(in)::t
       integer,intent(in):: nits
     end subroutine GS_Lex_const_Lap
     ! 
     subroutine Jacobi_const_Lap(phi,rhs,p,t,nits)
       use discretization
       use pms, only:verbose
       implicit none
       type(patcht)::p
       double precision,dimension(&
            p%all%lo%i:p%all%hi%i,&
            p%all%lo%j:p%all%hi%j,&
            p%all%lo%k:p%all%hi%k,nvar)::phi,rhs
       type(topot),intent(in)::t
       integer,intent(in):: nits
     end subroutine Jacobi_const_Lap
     !
     subroutine formExactU(exact,all,val,ip,dx)
       use discretization    
       implicit none
       type(box),intent(in)::val
       type(ipoint),intent(in)::ip
       type(box),intent(in)::all
       type(dpoint),intent(in):: dx
       double precision,intent(out)::exact(all%lo%i:all%hi%i,all%lo%j:all%hi%j,&
            all%lo%k:all%hi%k,1)
     end subroutine formExactU
     subroutine FormF(rhs,p,val,ip)
       use base_data_module 
       implicit none
       type(patcht)::p             ! max, dx & allocated
       type(box),intent(in)::val   ! needed to get correct coordinates with buffers in max (SR)
       type(ipoint),intent(in)::ip ! global position for coordinates
       double precision,intent(out)::rhs(p%all%lo%i:p%all%hi%i,p%all%lo%j:p%all%hi%j,&
            p%all%lo%k:p%all%hi%k,1)
     end subroutine FormF
     !
     double precision function norm(ux,all,val,dx,comm,type)
       use base_data_module
       use mpistuff
       implicit none 
       type(box),intent(in)::val
       type(box),intent(in)::all
       type(dpoint),intent(in):: dx
       double precision,intent(in)::ux(all%lo%i:all%hi%i,all%lo%j:all%hi%j,&
            all%lo%k:all%hi%k,1)
       integer,intent(in) :: type
       integer,intent(in)::comm
     end function norm
  end interface
contains
  !-----------------------------------------------------------------
  logical function is_top(g)
    use pms, only:bot_min_sz
    implicit none
    type(crs_patcht):: g
    is_top = ( &
         ((g%p%max%hi%i-g%p%max%lo%i+1).le.bot_min_sz) .or. &
         ((g%p%max%hi%j-g%p%max%lo%j+1).le.bot_min_sz) &
#ifndef TWO_D
         .or. ((g%p%max%hi%k-g%p%max%lo%k+1).le.bot_min_sz) &
#endif
         )
  end function is_top
  !-----------------------------------------------------------------
  subroutine new_grids(gc,gsr,NPeAxis,iPeAxis,nloc,comm3d)
    use pms,only:dom_max_sr_grids,dom_max_grids
    implicit none    
    type(ipoint)::nloc
    integer,intent(in):: comm3d
    integer :: NPeAxis(3),iPeAxis(3)
    type(crs_patcht),intent(out):: gc(0:dom_max_grids-1)
    type(sr_patcht),intent(out):: gsr(-dom_max_sr_grids:0)
    interface
       subroutine new_grids_private(gc,gsr,NPeAxis,iPeAxis,nloc,comm3d)
         use pms, only:dom_max_sr_grids,dom_max_grids
         use base_data_module
         implicit none
         type(ipoint)::nloc
         integer,intent(in)::comm3d
         integer::NPeAxis(3),iPeAxis(3)
         type(crs_patcht),intent(out):: gc(0:dom_max_grids-1)
         type(sr_patcht),intent(out):: gsr(-dom_max_sr_grids:0)
       end subroutine new_grids_private
    end interface
    call new_grids_private(gc,gsr,NPeAxis,iPeAxis,nloc,comm3d)
  end subroutine new_grids
  !-----------------------------------------------------------------
  subroutine destroy_grids(gc,gsr)
    use pms, only:dom_max_sr_grids,dom_max_grids
    use base_data_module
    implicit none    
    type(crs_patcht),intent(out):: gc(0:dom_max_grids-1)
    type(sr_patcht),intent(out):: gsr(-dom_max_sr_grids:0)
    interface       
       subroutine destroy_grids_private(gc,gsr)
         use pms, only:dom_max_sr_grids,dom_max_grids
         use base_data_module
         implicit none
         type(crs_patcht),intent(out):: gc(0:dom_max_grids-1)
         type(sr_patcht),intent(out):: gsr(-dom_max_sr_grids:0)
       end subroutine destroy_grids_private
    end interface
    call destroy_grids_private(gc,gsr)
  end subroutine destroy_grids
  !
end module grid_module
! 
module sr_interfaces
  use discretization
  interface
     subroutine copy_sr_crs2(crs1,crs2,sr1,sr2,cg0,srg0)
       use discretization
       implicit none
       type(sr_patcht),intent(in)::srg0
       type(crs_patcht),intent(in)::cg0
       double precision,intent(out),dimension(&
            cg0%p%all%lo%i:cg0%p%all%hi%i,& 
            cg0%p%all%lo%j:cg0%p%all%hi%j,& 
            cg0%p%all%lo%k:cg0%p%all%hi%k, nvar)::crs1,crs2
       double precision,intent(in),dimension(&
            srg0%p%all%lo%i:srg0%p%all%hi%i,&
            srg0%p%all%lo%j:srg0%p%all%hi%j,&
            srg0%p%all%lo%k:srg0%p%all%hi%k,nvar)::sr1,sr2
     end subroutine copy_sr_crs2
     !
     subroutine copy_crs_sr_w_bc_ex(u_sr,u_crs,srg0,cg0)
       use discretization
       implicit none
       type(sr_patcht),intent(in)::srg0
       type(crs_patcht),intent(in)::cg0
       double precision,intent(in),dimension(&
            cg0%p%all%lo%i:cg0%p%all%hi%i,& 
            cg0%p%all%lo%j:cg0%p%all%hi%j,& 
            cg0%p%all%lo%k:cg0%p%all%hi%k, nvar)::u_crs
       double precision,intent(out),dimension(&
            srg0%p%all%lo%i:srg0%p%all%hi%i,&
            srg0%p%all%lo%j:srg0%p%all%hi%j,&
            srg0%p%all%lo%k:srg0%p%all%hi%k,nvar)::u_sr
     end subroutine copy_crs_sr_w_bc_ex
  end interface
end module sr_interfaces

