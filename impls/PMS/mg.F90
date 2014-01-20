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
     subroutine SRFin(gsr,cg,fl,sr_ux,sr_rhs,sr_aux,Apply1,Apply2,Relax,Res0,errors)
       use pe_patch_data_module
       use pms, only:nvcycles,verbose,nsr,ncgrids,dom_max_sr_grids,dom_max_grids
       use error_data_module
       use mpistuff,only:mype
       implicit none
       type(sr_patcht),intent(in)::gsr(-dom_max_sr_grids:0)
       type(crs_patcht),intent(in)::cg(0:dom_max_grids-1)
       type(data_ptr),dimension(-nsr:0)::sr_ux,sr_rhs,sr_aux
       type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
       integer,intent(in)::fl
       double precision,intent(out)::Res0
       external Apply1, Apply2, Relax
     end subroutine SRFin
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
        sr_ux(0)%p => ux
        sr_rhs(0)%p => rhs
        sr_aux(0)%p => L2u_f
        do lev=-1,-nsr+1,-1
print *,'do mid SR for ',lev,lev+1
           NULLIFY(sr_ux(lev)%p,sr_ux(lev)%p,sr_aux(lev)%p)
           !call SRMid()
        end do
        lev = -nsr
print *,'do Final SR for',lev,lev+1        
        NULLIFY(sr_ux(lev)%p,sr_ux(lev)%p,sr_aux(lev)%p)
        call SRFin(gsr,cg,lev,sr_ux,sr_rhs,sr_aux,Apply1,Apply2,Relax,Res0,errors)
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

  ! output
  if (verbose.gt.10) then
     call formExactU(L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%t%ipe,cg(0)%p%dx)
     call WriteAVSFile(cg(0),ux,rhs,ux-L2u_f,1000)
  end if
  
  do lev=-nsr,-1,1
     DEALLOCATE(sr_ux(lev)%p,sr_rhs(lev)%p,sr_aux(lev)%p)
  end do
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
     if (verbose.gt.1) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') MGF bot |u_1|=',tt
     end if
  else
     ! go "down" the V, allating data (on stack), forming RHS, zero out u
     ! form RHS explicity on fine grid
     call formRHS(FC,g(lev+1)%p,g(lev+1)%p%max,g(lev+1)%t%ipe) 
     if (verbose.gt.2) then
        tt = norm(FC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'lev=',lev,') FMG: |f|=',tt
     end if
     call MGF(uxC,FC,g,lev+1,Apply1,Apply2,Relax,errors) 
     if(mype==0.and.verbose.gt.2)write(6,'(A,I2)') '[ 0] MGF: start lev=',lev
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
  
  cfoffset%i = 0;   cfoffset%j = 0;   cfoffset%k = 0; 

  if (verbose.gt.1) then  
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
       if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: after pre |u|=',tt
    end if
     !     restrict residual
     call Apply(Res,ux,g(lev)%p,g(lev)%t)
     Res = rhs - Res
    if (verbose.gt.2) then
       tt = norm(Res,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
       if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: after pre |r|=',tt
    end if
     !     restrict solution & RHS
     call RestrictFuse(g(lev+1),g(lev),g(lev+1)%t,cfoffset,FC,uxC,Res,ux,uxC)
     !     rhs = residual + Ac(Uc)
     if (verbose.gt.2) then
        tt = norm(FC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: after R |f|=',tt
        tt = norm(uxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: after R |u|=',tt
     end if
     call Apply(tmpC,uxC,g(lev+1)%p,g(lev+1)%t)
     FC = FC + tmpC
     if (verbose.gt.2) then
        tt = norm(FC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: coarse |r|=',tt
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
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: after MGV() |u|=',tt
     end if
     !     subtract off old Uc
     uxC = uxC - tmpC
     ! prolongate and correct
     call Prolong(ux,uxC,g(lev)%t,g(lev+1)%t,cfoffset,g(lev),g(lev+1))
     if (verbose.gt.2) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: after prol |u|=',tt
     end if
     ! end of v cycle, post smoothing
     call Relax(ux,rhs,g(lev)%p,g(lev)%t,nsmoothsup)
     if (verbose.gt.2) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if(mype==0)write(6,'(A,I2,A,E14.6)')'    lev=',lev,') VMG: after post smooth |u|=',tt
     end if
  end if
  return
end subroutine MGV

!-----------------------------------------------------------------------
! SR: SRFin - last leg of SR, collect functional, evanescent data
!-----------------------------------------------------------------------
subroutine SRFin(gsr,cg,lev,sr_ux,sr_rhs,sr_aux,Apply1,Apply2,Relax,Res0,errors)
  use pe_patch_data_module
  use pms
  use error_data_module
  use mpistuff,only:mype
  use discretization, only:nvar
  implicit none
  type(sr_patcht),intent(in)::gsr(-dom_max_sr_grids:0)
  type(crs_patcht),intent(in)::cg(0:dom_max_grids-1)
  type(data_ptr), dimension(-nsr:0)::sr_ux,sr_rhs,sr_aux
  type(error_data),intent(out)::errors(-nsr:ncgrids-1+nvcycles)
  integer,intent(in)::lev
  double precision,intent(out)::Res0
  external Apply1, Apply2, Relax
  ! 
  double precision,pointer,DIMENSION(:,:,:,:) :: uxF,rhsF,auxF,uxC,rhsC,auxC
  type(data_ptr),dimension(lev+1:0)::tmpC
  integer:: cl,iflv
  double precision:: norm,tt

  cl=lev+1
  print *,'SRFin: ASSOCIATED(uxF)=',ASSOCIATED(uxF)
  ALLOCATE(uxF(&
       gsr(lev)%p%all%lo%i:gsr(lev)%p%all%hi%i,&
       gsr(lev)%p%all%lo%j:gsr(lev)%p%all%hi%j,&
       gsr(lev)%p%all%lo%k:gsr(lev)%p%all%hi%k,nvar),rhsF(&
       gsr(lev)%p%all%lo%i:gsr(lev)%p%all%hi%i,&
       gsr(lev)%p%all%lo%j:gsr(lev)%p%all%hi%j,&
       gsr(lev)%p%all%lo%k:gsr(lev)%p%all%hi%k,nvar),auxF(&
       gsr(lev)%p%all%lo%i:gsr(lev)%p%all%hi%i,&
       gsr(lev)%p%all%lo%j:gsr(lev)%p%all%hi%j,&
       gsr(lev)%p%all%lo%k:gsr(lev)%p%all%hi%k,nvar))
  print *,'SRFin: rhsC=',rhsC(1,1,1,1)
  print *,'SRFin: fine grid:',&
       gsr(lev)%p%all%lo%i,gsr(lev)%p%all%hi%i,&
       gsr(lev)%p%all%lo%j,gsr(lev)%p%all%hi%j,&
       gsr(lev)%p%all%lo%k,gsr(lev)%p%all%hi%k,nvar
  ! 
  sr_ux(lev)%p  = uxF 
  sr_rhs(lev)%p = rhsF
  sr_aux(lev)%p = auxF
  uxC  =  sr_ux(cl)%p
  rhsC = sr_rhs(cl)%p
  auxC = sr_aux(cl)%p

  ! SR1 = Prol + SR2
  call formRHS(rhsF,gsr(lev)%p,gsr(lev)%p%max,gsr(lev)%t%ipe)
  uxF = 0.d0
  call Prolong(uxF,uxC,gsr(lev)%t,gsr(cl)%t,gsr(cl)%cfoffset,gsr(lev),gsr(cl))

  ! SR2 = smooth restrict
  print *,'SRFin: SR2 cl=',cl
  do iflv=lev,-1,1
     print *,'SRFin: SR2 with ',iflv,iflv+1     
     uxF  =  sr_ux(iflv)%p
     rhsF = sr_rhs(iflv)%p
     auxF = sr_aux(iflv)%p
     uxC  = sr_ux(iflv+1)%p
     rhsC = sr_rhs(iflv+1)%p
     auxC = sr_aux(iflv+1)%p
     ! SR2:    pre smoothing
     call Relax(uxF,rhsF,gsr(iflv)%p,gsr(iflv)%t,nsmoothsdown)
     !     restrict residual
     call Apply2(auxF,uxF,gsr(iflv)%p,gsr(iflv)%t)
     auxF = rhsF - auxF
     !     restrict solution & RHS
     call RestrictFuse(gsr(iflv+1),gsr(iflv),gsr(iflv+1)%t,gsr(iflv+1)%cfoffset,auxC,uxC,auxF,uxF,uxC)
     !     rhs = residual + Ac(Uc)
     ALLOCATE(tmpC(iflv+1)%p(&
          gsr(iflv+1)%p%all%lo%i:gsr(iflv+1)%p%all%hi%i,&
          gsr(iflv+1)%p%all%lo%j:gsr(iflv+1)%p%all%hi%j,&
          gsr(iflv+1)%p%all%lo%k:gsr(iflv+1)%p%all%hi%k,nvar))
     call Apply2(tmpC(iflv+1)%p,uxC,gsr(iflv+1)%p,gsr(iflv+1)%t)
     auxC = auxC + tmpC(iflv+1)%p
     !     Temporarily store uxC into tmpC
     tmpC(iflv+1)%p = uxC
  end do

  ! coarse grid solve
  call MGV(sr_ux(0)%p,sr_rhs(0)%p,cg,0,Apply1,Relax)

  ! SR3: prlongate, update, smooth 
  print *,'SRFin: SR3 cl=',cl
  do iflv=0,cl,-1
     print *,'SRFin: SR3 with ',iflv,iflv+1
     uxF  =  sr_ux(iflv)%p
     rhsF = sr_rhs(iflv)%p
     auxF = sr_aux(iflv)%p
     uxC = sr_ux(iflv+1)%p
     rhsC = sr_rhs(iflv+1)%p
     auxC = sr_aux(iflv+1)%p
     !     subtract off old Uc
     uxC = uxC - tmpC(iflv+1)%p
     DEALLOCATE(tmpC(iflv+1)%p)
     ! prolongate and correct
     call Prolong(uxF,uxC,gsr(iflv)%t,gsr(iflv+1)%t,gsr(iflv+1)%cfoffset,gsr(iflv),gsr(iflv+1))
     ! end of v cycle, post smoothing
     call Relax(uxF,rhsF,gsr(iflv)%p,gsr(iflv)%t,nsmoothsup)
  end do
  ! form errors & convergance measure, recursive so goes from coarse to fine (good)
  uxF  =  sr_ux(lev)%p
  rhsF = sr_rhs(lev)%p
  auxF = sr_aux(lev)%p
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

end subroutine SRFin
!----------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine WriteAVSFile(g,ux,rhs,error,index)
  use grid_module
  use pms
  use domain
  use mpistuff
  use iounits
  !use discretization
  implicit none
  type(crs_patcht),intent(in):: g
  !type(patcht)::p
  !type(box),intent(in)::val
  !type(ipoint),intent(in)::ip
  double precision,intent(in),dimension(&
       g%p%all%lo%i:g%p%all%hi%i,&
       g%p%all%lo%j:g%p%all%hi%j,&
       g%p%all%lo%k:g%p%all%hi%k,1)::ux,rhs,error
  integer,intent(in)::index
  !
  real,dimension(g%p%max%hi%i,g%p%max%hi%j,g%p%max%hi%k)::xn,yn,zn
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  integer:: nbytes,offset,itmp,nelements
  character*50 outfile,fldfile
  integer,parameter::ifld=91
  double precision::coord(3)

  ig = getIglobalx(g%p%max,g%t%ipe)
  jg = getIglobaly(g%p%max,g%t%ipe)
  kg = getIglobalz(g%p%max,g%t%ipe)
  do kk=1,g%p%max%hi%i
     do jj=1,g%p%max%hi%j
        do ii=1,g%p%max%hi%k
           xn(ii,jj,kk) = real(xl+(ig+ii-1)*g%p%dx%i-0.5*g%p%dx%i)
           yn(ii,jj,kk) = real(yl+(jg+jj-1)*g%p%dx%j-0.5*g%p%dx%j)
#ifndef TWO_D
           zn(ii,jj,kk) = real(zl+(kg+kk-1)*g%p%dx%k-0.5*g%p%dx%k)
#endif
        end do
     end do
  end do

  !     File name for data 
  write(outfile,1000) float(mype)/1000.0 !,float(index)/1000000.0
1000 format('u_rhs_error',f4.3,'.dat') !,f7.6)
  !     Write out fld file
  if(g%t%ipe%i.eq.1.and.g%t%ipe%j.eq.1.and.g%t%ipe%k.eq.1) then
     write(fldfile,2000) float(index)/1000000.0
2000 format('fbov',f7.6,'.bov')
     open(ifld,file=fldfile,form='formatted')
     write(ifld,*) 'TIME: 0.0'
     write(ifld,*) 'DATA_FILE: u_rhs_error.%3d.dat'
     write(ifld,*) 'DATA_SIZE: ',g%p%max%hi%i*g%t%npe%i+1,g%p%max%hi%j*g%t%npe%j+1,g%p%max%hi%k*g%t%npe%k+1
     write(ifld,*) 'DATA_FORMAT: DOUBLE'
     write(ifld,*) 'VARIABLE: u,rhs,error'
     !     How to detect little vs. big endian automatically in Fortran?
     write(ifld,*) 'DATA_ENDIAN: LITTLE'
     write(ifld,*) 'CENTERING: zonal'
     write(ifld,*) 'BYTE_OFFSET: 4'
     write(ifld,*) 'BRICK_ORIGIN: 0.0 0.0 0.0'
     write(ifld,*) 'BRICK_SIZE: 1.0 1.0 1.0'
     close (ifld)
  endif

  open(itecoutput,file=outfile,form='unformatted')  
  do ii=1,g%p%max%hi%i
     do jj=1,g%p%max%hi%j
        do kk=1,g%p%max%hi%k
           write(itecoutput) ux(ii,jj,kk,1)
           write(itecoutput) rhs(ii,jj,kk,1)
           write(itecoutput) error(ii,jj,kk,1)
        end do
     end do
  end do
  close(itecoutput)

end subroutine WriteAVSFile
