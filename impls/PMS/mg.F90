!-----------------------------------------------------------------
!     Ravi Samtaney & Mark Adams
!     Copyright 2014
!-----------------------------------------------------------------
subroutine solve(cg,gsr,Apply1,Apply2,Relax,Res0,errors,nViters)
  use pe_patch_data_module
  use pms, only:dom_max_grids,nvcycles,nfcycles,rtol,verbose,dom_max_sr_grids,nsr,ncgrids
  use mpistuff
  use error_data_module
  use discretization, only:nvar
  implicit none
  !
  type(crs_patcht),intent(in):: cg(0:dom_max_grids-1)
  type(sr_patcht),intent(in):: gsr(-dom_max_sr_grids:0)
  type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
  external Apply1, Apply2, Relax
  double precision,intent(out)::Res0
  integer,intent(out):: nViters
  !     Local Vars
  integer:: iter,lev
  double precision:: Res,Reslast
  double precision,target,dimension(& ! alloc fine grid data here
       cg(0)%p%all%lo%i:cg(0)%p%all%hi%i,&
       cg(0)%p%all%lo%j:cg(0)%p%all%hi%j,&
       cg(0)%p%all%lo%k:cg(0)%p%all%hi%k,nvar)::L2u_f,rhs,ux
  double precision norm
  type(data_ptr),dimension(-nsr:0)::sr_ux,sr_rhs,sr_aux
  interface
     subroutine MGSR(gsr,cg,fl,sr_ux,sr_rhs,sr_aux,uxC0,auxC0,Apply1,Apply2,Relax,Res0,errors)
       use pe_patch_data_module
       use pms, only:nvcycles,verbose,nsr,ncgrids,dom_max_sr_grids,dom_max_grids
       use error_data_module
       use mpistuff,only:mype
       use discretization, only:nvar
       implicit none
       type(sr_patcht),intent(in)::gsr(-dom_max_sr_grids:0)
       type(crs_patcht),intent(in)::cg(0:dom_max_grids-1)
       type(data_ptr),dimension(-nsr:0)::sr_ux,sr_rhs,sr_aux
       double precision,intent(in),dimension(&
            cg(0)%p%all%lo%i:cg(0)%p%all%hi%i,& 
            cg(0)%p%all%lo%j:cg(0)%p%all%hi%j,& 
            cg(0)%p%all%lo%k:cg(0)%p%all%hi%k,nvar)::auxC0,uxC0
       type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
       integer,intent(in)::fl
       double precision,intent(out)::Res0
       external Apply1, Apply2, Relax
     end subroutine MGSR
  end interface

  call formRHS(rhs,cg(0)%p,cg(0)%p%max,cg(0)%t%ipe)
  Res0 = norm(rhs,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
  if(verbose.gt.1) then
     call formExactU(L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%t%ipe,cg(0)%p%dx)
     Reslast = norm(L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2) 
     if(mype==0) write(6,'(A,E14.6,A,E14.6,A,I5)') '0) solve: |f|_2=',&
          Res0,', |u_exact|_2=',Reslast,', nx=',(cg(0)%p%max%hi%i-cg(0)%p%max%lo%i+1)*cg(0)%t%npe%i
  end if
  Reslast = Res0
  ! FMG (SR coarse grid) solve
  iter = 0
  err_lev = 0
  if (nfcycles .ne. 0) then
     iter = iter + 1
     ! SR start of FMG
     call MGF(ux,rhs,cg,0,Apply1,Apply2,Relax,errors)
     ! diagnostics
     call Apply2(L2u_f,ux,cg(0)%p,cg(0)%t)
     Res = norm(rhs-L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
     if(mype==0.and.verbose.gt.0)write(6,'(I2,A,E14.6,A,F10.6)') &
          iter,') solve: FMG done |f-Au|_2=',Res,', rate=',Res/Reslast
     Reslast=Res
     ! SR 
     if (nsr.gt.0) then
        NULLIFY(sr_ux(0)%ptr,sr_ux(0)%ptr,sr_aux(0)%ptr)
        do lev=-1,-nsr,-1
           NULLIFY(sr_ux(lev)%ptr,sr_ux(lev)%ptr,sr_aux(lev)%ptr)
           call MGSR(gsr,cg,lev,sr_ux,sr_rhs,sr_aux,ux,L2u_f,Apply1,Apply2,Relax,Res0,errors)
        end do
     end if
  else
     ! initailize what is done in FMG
     call formRHS(rhs,cg(0)%p,cg(0)%p%max,cg(0)%t%ipe) ! form RHS explicity on fine grid for vycycles
     Res = Res0
     ux=0.d0
     if (nsr.ne.0) stop 'SR but no F-cycle'
  end if
 
  ! finish with V-cycle - no V with SR
  nViters = 0
  do while(Res/Res0 > rtol .and. iter < nvcycles .and. nsr==0)
     call MGV(ux,rhs,cg,0,Apply1,Relax)
     iter=iter+1 ! 2 for FMG and 1 for V 
     ! residual
     call Apply2(L2u_f,ux,cg(0)%p,cg(0)%t)            
     Res = norm(rhs-L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
     if(mype==0.and.verbose.gt.0)write(6,'(I2,A,E14.6,A,F10.6)') &
          iter,') solve: |f-Au|_2=',Res,', rate=',Res/Reslast     
     ! form errors & convergance measure
     errors(err_lev)%resid = Res
     call formExactU(L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%t%ipe,cg(0)%p%dx)
     errors(err_lev)%uerror = norm(L2u_f-ux,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,error_norm)
     err_lev = err_lev + 1
     Reslast=Res
     nViters = nViters + 1
  enddo
  
  if (nsr.gt.0) then
     do lev=-nsr,0,1
        DEALLOCATE(sr_ux(lev)%ptr,sr_rhs(lev)%ptr,sr_aux(lev)%ptr)
     end do
  end if
  return
end subroutine solve
!-----------------------------------------------------------------------
! MGF: recursive F-cycle with data on stack
!-----------------------------------------------------------------------
recursive subroutine MGF(ux,rhs,g,lev,Apply1,Apply2,Relax,errors)
  use grid_module
  use mpistuff
  use error_data_module
  use pms
  use discretization, only:nvar
  implicit none
  integer,intent(in) :: lev
  type(crs_patcht),intent(in):: g(0:dom_max_grids-1)
  type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
  double precision,intent(in),dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,& 
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,& 
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k, nvar) :: rhs
  double precision,intent(out),dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,&
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,& 
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k, nvar) :: ux
  external Apply1, Apply2, Relax
  ! Local Vars
  integer :: ii,jj,kk,ist,jst,kst
  double precision :: norm,tt
  double precision,dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,&
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,&
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k,nvar)::tmpF
  ! Coarse data. This is new data on the stack
  double precision,dimension(&
       g(lev+1)%p%all%lo%i:g(lev+1)%p%all%hi%i,&
       g(lev+1)%p%all%lo%j:g(lev+1)%p%all%hi%j,&
       g(lev+1)%p%all%lo%k:g(lev+1)%p%all%hi%k, nvar):: uxC,FC
  type(ipoint)::cfoffset

  cfoffset%i = 0;   cfoffset%j = 0;   cfoffset%k = 0; 

  ux=0.d0 ! everthing uses initial solution except FMG
  if( is_top(g(lev)) ) then
     ! the real start of FMG - coarse grid solve
     if (lev.gt.0 .and. is_top(g(lev-1))) stop 'MGF: not coarsest grid'
     if(mype==0.and.verbose.gt.2)write(6,'(A,I2)') '[ 0]MGF: bot solve lev=',lev
     call Relax(ux,rhs,g(lev)%p,g(lev)%t,ncoarsesolveits)  ! coar grid solve
     if (verbose.gt.3) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') MGF bot |u_1|=',tt
     end if
  else
     ! go "down" the V, allating data (on stack), forming RHS, zero out u
     ! form RHS explicity on fine grid
     call formRHS(FC,g(lev+1)%p,g(lev+1)%p%max,g(lev+1)%t%ipe) 
     if (verbose.gt.3) then
        tt = norm(FC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') F: |f|=',tt
     end if
     call MGF(uxC,FC,g,lev+1,Apply1,Apply2,Relax,errors) 
     if(mype==0.and.verbose.gt.2)write(6,'(A,I2)') '[ 0] M: start lev=',lev
     ! prolongate and correct
     call Prolong(ux,uxC,g(lev)%t,g(lev+1)%t,cfoffset,g(lev),g(lev+1))
     if (verbose.gt.2) then
        tt = norm(uxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev+1,') MGF prol src |u|=',tt
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') MGF prol dest |u|=',tt
     end if

     ! start of multiple v-cycle
     jj = nfmgvcycles !; if(lev==0) jj=10 ! accurate coarse solve
     do ii=1,jj
        call MGV(ux,rhs,g,lev,Apply1,Relax)
     enddo
  end if

  ! form errors & convergance measure, recursive so goes from coarse to fine (good)
  call formExactU(tmpF,g(lev)%p%all,g(lev)%p%max,g(lev)%t%ipe,g(lev)%p%dx)
  errors(err_lev)%uerror = norm(tmpF-ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,error_norm)
  call Apply2(tmpF,ux,g(lev)%p,g(lev)%t)
  errors(err_lev)%resid = norm(tmpF-rhs,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)

  if (verbose.gt.0) then
     if(mype==0)write(6,'(A,I2,A,E14.6,A,E14.6,A,I8,I8,I8,A,I4,I4,I4)') &
          '     lev=',lev,') MGF |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror,&
          ', n=',(g(lev)%p%max%hi%i-g(lev)%p%max%lo%i+1)*g(lev)%t%npe%i,&
                 (g(lev)%p%max%hi%j-g(lev)%p%max%lo%j+1)*g(lev)%t%npe%j,&
                 (g(lev)%p%max%hi%k-g(lev)%p%max%lo%k+1)*g(lev)%t%npe%k,&
          ', npe=',g(lev)%t%npe%i,g(lev)%t%npe%j,g(lev)%t%npe%k
  end if
  
  err_lev = err_lev + 1
  return
end subroutine MGF
!-----------------------------------------------------------------------
recursive subroutine MGV(ux,rhs,g,lev,Apply,Relax)
  use grid_module
  use mpistuff
  use pms
  use discretization
  implicit none
  integer,intent(in)::lev
  type(crs_patcht),intent(in)::g(0:dom_max_grids-1)
  double precision,dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i, g(lev)%p%all%lo%j:g(lev)%p%all%hi%j, &
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k,nvar) :: ux, rhs
  external Apply,Relax
  !
  type(ipoint)::cfoffset
  double precision norm,tt
  integer:: ii,jj
  !     Coarse, here is the allocation on stack
  double precision,dimension(&
       g(lev+1)%p%all%lo%i:g(lev+1)%p%all%hi%i,&
       g(lev+1)%p%all%lo%j:g(lev+1)%p%all%hi%j, &
       g(lev+1)%p%all%lo%k:g(lev+1)%p%all%hi%k,nvar) :: uxC, tmpC, FC, ResC
  ! This level
  double precision,dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,&
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,& 
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k, nvar) :: Res
  
  ! non-SR need a dummy for the offset
  cfoffset%i = 0;   cfoffset%j = 0;   cfoffset%k = 0; 

  if (verbose.gt.2) then  
     tt = norm( ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: u_0=',tt
     tt = norm(rhs,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6)')'       lev=',lev,') V: f_0=',tt
  end if
  if( is_top(g(lev)) ) then
     call Relax(ux,rhs,g(lev)%p,g(lev)%t,ncoarsesolveits)
  else
     !     pre smoothing
     call Relax(ux,rhs,g(lev)%p,g(lev)%t,nsmoothsdown)
     if (verbose.gt.2) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: after pre |u|=',tt
     end if
     !     residual
     call Apply(Res,ux,g(lev)%p,g(lev)%t)
     Res = rhs - Res
     if (verbose.gt.2) then
        tt = norm(Res,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: after pre |r|=',tt
     end if
     !     restrict solution & RHS
     call RestrictFuse(g(lev+1),g(lev),g(lev+1)%t,cfoffset,FC,uxC,Res,ux)
     !     rhs = residual + Ac(Uc)
     if (verbose.gt.2) then
        tt = norm(FC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: after R |f|=',tt
        tt = norm(uxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: after R |u|=',tt
     end if
     call Apply(tmpC,uxC,g(lev+1)%p,g(lev+1)%t)
     FC = FC + tmpC
     if (verbose.gt.2) then
        tt = norm(FC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: coarse |r|=',tt
     end if
     !     Temporarily store uxC into tmpC
     tmpC = uxC
     !     start of v(1) or w(2) cycle 
     jj = ncycles ! ; if(lev==0) jj=40 ! two level method, sort of
     do ii=1,jj
        call MGV(uxC,FC,g,lev+1,Apply,Relax)
     enddo
     if (verbose.gt.2) then
        tt = norm(uxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: after MGV() |u|=',tt
     end if
     !     subtract off old Uc
     tmpC = uxC - tmpC
     ! prolongate and correct
     call Prolong(ux,tmpC,g(lev)%t,g(lev+1)%t,cfoffset,g(lev),g(lev+1))
     if (verbose.gt.2) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: after prol |u|=',tt
     end if
     ! end of v cycle, post smoothing
     call Relax(ux,rhs,g(lev)%p,g(lev)%t,nsmoothsup)
     if (verbose.gt.2) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') V: after post smooth |u|=',tt
     end if
  end if
  return
end subroutine MGV
!-----------------------------------------------------------------------
! SR: copy_sr_crs2 -- copy two vectors from SR vectors (sr0) to normal (cg0)
!-----------------------------------------------------------------------
subroutine copy_sr_crs2(crs1,crs2,sr1,sr2,cg0,gsr0)
  use pe_patch_data_module
  use pms
  use mpistuff,only:mype
  use discretization, only:nvar
  implicit none
  type(sr_patcht),intent(in)::gsr0
  type(crs_patcht),intent(in)::cg0
  double precision,intent(out),dimension(&
       cg0%p%all%lo%i:cg0%p%all%hi%i,& 
       cg0%p%all%lo%j:cg0%p%all%hi%j,& 
       cg0%p%all%lo%k:cg0%p%all%hi%k, nvar)::crs1,crs2
  double precision,intent(in),dimension(&
       gsr0%p%all%lo%i:gsr0%p%all%hi%i,&
       gsr0%p%all%lo%j:gsr0%p%all%hi%j,&
       gsr0%p%all%lo%k:gsr0%p%all%hi%k,nvar)::sr1,sr2
  !
  integer::ii,jj,kk,sri,srj,srk
  ! copy
  ii=0
  do sri=gsr0%val%lo%i,gsr0%val%hi%i
     ii=ii+1
     jj=0
     do srj=gsr0%val%lo%j,gsr0%val%hi%j
        jj=jj+1
        kk=0
        do srk=gsr0%val%lo%k,gsr0%val%hi%k
           kk=kk+1
           crs1(ii,jj,kk,:) = sr1(sri,srj,srk,:)
           crs2(ii,jj,kk,:) = sr2(sri,srj,srk,:)
        end do
     end do
  end do
  print *,'copy_sr_crs2: id should be equal i=',ii,cg0%p%max%hi%i,', j=',jj,cg0%p%max%hi%j,', k=',kk,cg0%p%max%hi%k
  return
end subroutine copy_sr_crs2
!-----------------------------------------------------------------------
! SR: copy_crs_sr_w_bc_ex - copy vector from normal to SR
!-----------------------------------------------------------------------
subroutine copy_crs_sr_w_bc_ex(u_sr,u_crs,gsr0,cg0)
  use pe_patch_data_module
  use pms
  use mpistuff,only:mype
  use discretization
  implicit none
  type(sr_patcht),intent(in)::gsr0
  type(crs_patcht),intent(in)::cg0
  double precision,intent(in),dimension(&
       cg0%p%all%lo%i:cg0%p%all%hi%i,& 
       cg0%p%all%lo%j:cg0%p%all%hi%j,& 
       cg0%p%all%lo%k:cg0%p%all%hi%k, nvar)::u_crs
  double precision,intent(out),dimension(&
       gsr0%p%all%lo%i:gsr0%p%all%hi%i,&
       gsr0%p%all%lo%j:gsr0%p%all%hi%j,&
       gsr0%p%all%lo%k:gsr0%p%all%hi%k,nvar)::u_sr
  !
  integer::ii,jj,kk,sri,srj,srk,nsg_orig,nbuff
  type(patcht)::p
  type(ipoint)::shift
  ! copy to valid from whole
  ii=0
  do sri=gsr0%val%lo%i,gsr0%val%hi%i
     ii=ii+1
     jj=0
     do srj=gsr0%val%lo%j,gsr0%val%hi%j
        jj=jj+1
        kk=0
        do srk=gsr0%val%lo%k,gsr0%val%hi%k
           kk=kk+1
           u_sr(sri,srj,srk,:) = u_crs(ii,jj,kk,:)
        end do
     end do
  end do

  ! we have to hack into this to get the exchange to work
  nsg_orig = nsg%i
  if(nsg%j.ne.nsg_orig.or.nsg%k.ne.nsg_orig) stop 'copy_crs_sr_w_bc_ex: nsg%j.ne.nsg_orig'
  ! set fake ghost sizes
  if (gsr0%p%all%hi%i - gsr0%p%max%hi%i .gt. nsg%i) nsg%i = gsr0%p%all%hi%i - gsr0%p%max%hi%i
  if (-gsr0%p%all%hi%i+1.gt.nsg%i) nsg%i = -gsr0%p%all%hi%i+1
  ! 
  if (gsr0%p%all%hi%j - gsr0%p%max%hi%j .gt. nsg%j) nsg%j = gsr0%p%all%hi%j - gsr0%p%max%hi%j
  if (-gsr0%p%all%hi%j+1.gt.nsg%j) nsg%j = -gsr0%p%all%hi%j+1
  ! 
  if (gsr0%p%all%hi%k - gsr0%p%max%hi%k .gt. nsg%k) nsg%k = gsr0%p%all%hi%k - gsr0%p%max%hi%k
  if (-gsr0%p%all%hi%k+1.gt.nsg%k) nsg%k = -gsr0%p%all%hi%k+1

  ! strange patch with big exchange areas
  p%dx%i = cg0%p%dx%i ! just copy dx
  p%dx%j = cg0%p%dx%j
  p%dx%k = cg0%p%dx%k
  ! shift alloc
  shift%i=gsr0%val%lo%i-1
  shift%j=gsr0%val%lo%j-1
  shift%k=gsr0%val%lo%k-1
  
  ! shift alloc region to put valid in max
  p%all%lo%i=gsr0%p%all%lo%i-shift%i
  p%all%hi%i=gsr0%p%all%hi%i-shift%i
  p%all%lo%j=gsr0%p%all%lo%j-shift%j
  p%all%hi%j=gsr0%p%all%hi%j-shift%j
#ifndef TWO_D
  p%all%lo%k=gsr0%p%all%lo%k-shift%k
  p%all%hi%k=gsr0%p%all%hi%k-shift%k
#else
  p%all%lo%k=1
  p%all%hi%k=1
#endif 
  ! max is valid size
  p%max%hi%i=gsr0%val%hi%i-gsr0%val%lo%i+1
  p%max%hi%j=gsr0%val%hi%j-gsr0%val%lo%j+1
#ifndef TWO_D
  p%max%hi%k=gsr0%val%hi%k-gsr0%val%lo%k+1
#else
  p%max%hi%k=1
#endif 

  call SetBCs(u_sr,p,gsr0%t)

  nsg%i = nsg_orig
  nsg%j = nsg_orig
  nsg%k = nsg_orig
  return
end subroutine copy_crs_sr_w_bc_ex
!-----------------------------------------------------------------------
! SR: MGSR - last leg of SR, collect functional, evanescent data
!-----------------------------------------------------------------------
subroutine MGSR(gsr,cg,lev,sr_ux,sr_rhs,sr_aux,uxC0,auxC0,Apply1,Apply2,Relax,Res0,errors)
  use pe_patch_data_module
  use pms
  use error_data_module
  use mpistuff,only:mype
  use discretization, only:nvar
  implicit none
  type(sr_patcht),intent(in)::gsr(-dom_max_sr_grids:0)
  type(crs_patcht),intent(in)::cg(0:dom_max_grids-1)
  type(data_ptr),dimension(-nsr:0)::sr_ux,sr_rhs,sr_aux
  double precision,dimension(&
       cg(0)%p%all%lo%i:cg(0)%p%all%hi%i,& 
       cg(0)%p%all%lo%j:cg(0)%p%all%hi%j,& 
       cg(0)%p%all%lo%k:cg(0)%p%all%hi%k,nvar)::uxC0,auxC0 ! in/out, buffer
!!$  double precision,dimension(&
!!$       cg(0)%p%all%lo%i:cg(0)%p%all%hi%i,&
!!$       cg(0)%p%all%lo%j:cg(0)%p%all%hi%j,& 
!!$       cg(0)%p%all%lo%k:cg(0)%p%all%hi%k,nvar)::uxC0,resC0 ! buffers for coarse grid V-cycle
  type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
  integer,intent(in)::lev ! negative index into gsr(lev)
  double precision,intent(out)::Res0
  external Apply1, Apply2, Relax
  ! 
  double precision,pointer,DIMENSION(:,:,:,:)::uxF,rhsF,auxF,uxC,auxC
  type(data_ptr),dimension(lev+1:0)::tmpSRC
  integer::iflv,AllocateStatus
  double precision:: norm,tt,t2

  ! allocate coarse grid if needed, first call
  if (.not.ASSOCIATED(sr_ux(0)%ptr)) then
     if (lev.ne.-1) stop 'MGSR: lev.ne.-1'
     ALLOCATE(sr_ux(0)%ptr(&
       gsr(0)%p%all%lo%i:gsr(0)%p%all%hi%i,&
       gsr(0)%p%all%lo%j:gsr(0)%p%all%hi%j,&
       gsr(0)%p%all%lo%k:gsr(0)%p%all%hi%k,nvar),sr_rhs(0)%ptr(&
       gsr(0)%p%all%lo%i:gsr(0)%p%all%hi%i,&
       gsr(0)%p%all%lo%j:gsr(0)%p%all%hi%j,&
       gsr(0)%p%all%lo%k:gsr(0)%p%all%hi%k,nvar),sr_aux(0)%ptr(&
       gsr(0)%p%all%lo%i:gsr(0)%p%all%hi%i,&
       gsr(0)%p%all%lo%j:gsr(0)%p%all%hi%j,&
       gsr(0)%p%all%lo%k:gsr(0)%p%all%hi%k,nvar), STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
     ! copy in solution, apply BC, and exchange
     uxC => sr_ux(0)%ptr 
     call copy_crs_sr_w_bc_ex(uxC,uxC0,gsr(0),cg(0)) ! input


     tt = norm(uxC0,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') MGSR copy_crs_sr_w_bc_ex |u_src|=',tt
     t2 = norm(uxC,gsr(0)%p%all,gsr(0)%val,gsr(0)%p%dx,gsr(0)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6,E14.6)')'lev=',lev,') MGSR copy_crs_sr_w_bc_ex |u_dest|=',t2,tt-t2
     if(abs(tt-t2)/t2>1.d-12) stop 'MGSR: big error'


     ! just compute RHS
     uxC => sr_rhs(0)%ptr
     call formRHS(uxC,gsr(0)%p,gsr(0)%p%max,gsr(0)%t%ipe)
  end if

  ! allocate new fine grid
  ALLOCATE(sr_ux(lev)%ptr(&
       gsr(lev)%p%all%lo%i:gsr(lev)%p%all%hi%i,&
       gsr(lev)%p%all%lo%j:gsr(lev)%p%all%hi%j,&
       gsr(lev)%p%all%lo%k:gsr(lev)%p%all%hi%k,nvar),sr_rhs(lev)%ptr(&
       gsr(lev)%p%all%lo%i:gsr(lev)%p%all%hi%i,&
       gsr(lev)%p%all%lo%j:gsr(lev)%p%all%hi%j,&
       gsr(lev)%p%all%lo%k:gsr(lev)%p%all%hi%k,nvar),sr_aux(lev)%ptr(&
       gsr(lev)%p%all%lo%i:gsr(lev)%p%all%hi%i,&
       gsr(lev)%p%all%lo%j:gsr(lev)%p%all%hi%j,&
       gsr(lev)%p%all%lo%k:gsr(lev)%p%all%hi%k,nvar),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  ! SR1 = Prol + SR2
  rhsF => sr_rhs(lev)%ptr
  call formRHS(rhsF,gsr(lev)%p,gsr(lev)%p%max,gsr(lev)%t%ipe)
  uxC => sr_ux(lev+1)%ptr 
  uxF => sr_ux(lev)%ptr 
  uxF = 0.d0
  call Prolong(uxF,uxC,gsr(lev)%t,gsr(lev+1)%t,gsr(lev+1)%cfoffset,gsr(lev),gsr(lev+1))
  if (verbose.gt.2) then
     tt = norm(uxC,gsr(lev+1)%p%all,gsr(lev+1)%val,gsr(lev+1)%p%dx,gsr(lev+1)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev+1,') SR FMG prol src |u|=',tt
     tt = norm(uxF,gsr(lev)%p%all,gsr(lev)%val,gsr(lev)%p%dx,gsr(lev)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') SR FMG prol dest |u|=',tt
  end if
  
  ! SR2 = smooth restrict
  do iflv=lev,-1,1
     uxF  => sr_ux (iflv)%ptr
     rhsF => sr_rhs(iflv)%ptr
     auxF => sr_aux(iflv)%ptr
     uxC  => sr_ux (iflv+1)%ptr 
     !rhsC => sr_rhs(iflv+1)%ptr
     auxC => sr_aux(iflv+1)%ptr
     ! SR2:    pre smoothing
print *,'MGSR: call relax ',iflv
     call Relax(uxF,rhsF,gsr(iflv)%p,gsr(iflv)%t,nsmoothsdown)
print *,'MGSR: call relax DONE',iflv
     if (verbose.gt.2) then
        tt = norm(uxF,gsr(iflv)%p%all,gsr(iflv)%val,gsr(iflv)%p%dx,gsr(iflv)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after pre |u|=',tt
     end if
     !     residual
     call Apply2(auxF,uxF,gsr(iflv)%p,gsr(iflv)%t)
     auxF = rhsF - auxF
     if (verbose.gt.2) then
        tt = norm(auxF,gsr(iflv)%p%all,gsr(iflv)%val,gsr(iflv)%p%dx,gsr(iflv)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after pre |r|=',tt
     end if
     !     restrict solution & RHS
     call RestrictFuse(gsr(iflv+1),gsr(iflv),gsr(iflv+1)%t,gsr(iflv+1)%cfoffset,auxC,uxC,auxF,uxF)
     if (verbose.gt.2) then
        tt = norm(auxF,gsr(iflv+1)%p%all,gsr(iflv+1)%val,gsr(iflv+1)%p%dx,gsr(iflv+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after R |f|=',tt
        tt = norm(uxC,gsr(iflv+1)%p%all,gsr(iflv+1)%val,gsr(iflv+1)%p%dx,gsr(iflv+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after R |u|=',tt
     end if
     !     rhs = residual + Ac(Uc)
     ALLOCATE(tmpSRC(iflv+1)%ptr(&
          gsr(iflv+1)%p%all%lo%i:gsr(iflv+1)%p%all%hi%i,&
          gsr(iflv+1)%p%all%lo%j:gsr(iflv+1)%p%all%hi%j,&
          gsr(iflv+1)%p%all%lo%k:gsr(iflv+1)%p%all%hi%k,nvar),STAT=AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
     call Apply2(tmpSRC(iflv+1)%ptr,uxC,gsr(iflv+1)%p,gsr(iflv+1)%t)
     auxC = auxC + tmpSRC(iflv+1)%ptr ! FC in MGV
     if (verbose.gt.2) then
        tt = norm(auxC,gsr(iflv+1)%p%all,gsr(iflv+1)%val,gsr(iflv+1)%p%dx,gsr(iflv+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: coarse |r|=',tt
     end if
     !     Temporarily store uxC into tmpSRC
     tmpSRC(iflv+1)%ptr = uxC ! level 0 at end
  end do

  ! copy U & F from
  print *,'1) auxC(1111)',auxC(1,1,1,1)
  auxC => sr_aux(0)%ptr ! FC in MGV
  print *,'2) auxC(1111)',auxC(1,1,1,1)
  print *,'1) uxC(1111)',uxC(1,1,1,1)
  uxC => tmpSRC(0)%ptr
  print *,'2) uxC(1111)',uxC(1,1,1,1)
  

  tt = norm(auxC,gsr(0)%p%all,gsr(0)%val,gsr(0)%p%dx,gsr(0)%t%comm,2)
  if(mype==0)write(6,'(A,I2,A,E14.6)')'***1    lev=',lev,') SR: 0 |r|=',tt
  

  call copy_sr_crs2(uxC0,auxC0,uxC,auxC,cg(0),gsr(0))

  tt = norm(auxC0,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
  if(mype==0)write(6,'(A,I2,A,E14.6)')'***2    lev=',lev,') SR: 0 |r|=',tt
  

  ! coarse grid solve - these pointers are coarse data types, no SR buffers
  call MGV(uxC0,auxC0,cg,0,Apply1,Relax)
  ! copy coarse grid do SR coarse grid (cg.0-->sr.0)  

  if (verbose.gt.2) then
     tt = norm(uxC,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after MGV() |u_0_crs|=',tt
  end if

  uxC => sr_ux(0)%ptr 
  call copy_crs_sr_w_bc_ex(uxC,uxC0,gsr(0),cg(0))

  if (verbose.gt.2) then
     tt = norm(uxC,gsr(0)%p%all,gsr(0)%val,gsr(0)%p%dx,gsr(0)%t%comm,2)
     if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after MGV() |u_0_sr|=',tt
  end if


  ! SR3: prlongate, update, smooth 
  print *,'MGSR: SR3 cl=',lev+1
  do iflv=-1,lev,-1
     print *,'MGSR: SR3 with ',iflv,iflv+1
     uxF  =>  sr_ux(iflv)%ptr
     rhsF => sr_rhs(iflv)%ptr
     auxF => sr_aux(iflv)%ptr
     uxC => sr_ux(iflv+1)%ptr
     !rhsC => sr_rhs(iflv+1)%ptr
     auxC => sr_aux(iflv+1)%ptr
     !     subtract off old Uc
     uxC = uxC - tmpSRC(iflv+1)%ptr ! increment now
     DEALLOCATE(tmpSRC(iflv+1)%ptr)
     ! prolongate and correct
     call Prolong(uxF,uxC,gsr(iflv)%t,gsr(iflv+1)%t,gsr(iflv+1)%cfoffset,gsr(iflv),gsr(iflv+1))
     if (verbose.gt.2) then
        tt = norm(uxF,gsr(iflv)%p%all,gsr(iflv)%val,gsr(iflv)%p%dx,gsr(iflv)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after prol |u|=',tt
     end if
     ! end of v cycle, post smoothing
     call Relax(uxF,rhsF,gsr(iflv)%p,gsr(iflv)%t,nsmoothsup)
     if (verbose.gt.2) then
        tt = norm(uxF,gsr(iflv)%p%all,gsr(iflv)%val,gsr(iflv)%p%dx,gsr(iflv)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') SR: after post smooth |u|=',tt
     end if
  end do

  ! form errors & convergance measure, recursive so goes from coarse to fine (good)
  call formExactU(auxF,gsr(lev)%p%all,gsr(lev)%val,gsr(lev)%t%ipe,gsr(lev)%p%dx)
  errors(err_lev)%uerror=norm(auxF-uxF,gsr(lev)%p%all,gsr(lev)%val,gsr(lev)%p%dx,gsr(lev)%t%comm,error_norm)
  call Apply2(auxF,uxF,gsr(lev)%p,gsr(lev)%t)
  errors(err_lev)%resid= norm(auxF-uxF,gsr(lev)%p%all,gsr(lev)%val,gsr(lev)%p%dx,gsr(lev)%t%comm,2)

  if (verbose.gt.0) then
     if(mype==0)write(6,'(A,I2,A,E14.6,A,E14.6,A,I8,I8,I8,A,I4,I4,I4)') &
          '     lev=',lev,') SR |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror,', n loc =',&
          (gsr(lev)%p%max%hi%i-gsr(lev)%p%max%lo%i+1),&
          (gsr(lev)%p%max%hi%j-gsr(lev)%p%max%lo%j+1),&
          (gsr(lev)%p%max%hi%k-gsr(lev)%p%max%lo%k+1)
  end if  
  err_lev = err_lev + 1

end subroutine MGSR
