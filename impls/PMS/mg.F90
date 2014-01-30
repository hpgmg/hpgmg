!-----------------------------------------------------------------
!     Ravi Samtaney & Mark Adams
!     Copyright 2014
!-----------------------------------------------------------------
subroutine solve(cg,srg,Apply1,Apply2,Relax,Res0,errors,nViters)
  use grid_module
  use pms, only:dom_max_grids,nvcycles,nfcycles,rtol,verbose,dom_max_sr_grids,nsr,ncgrids
  use mpistuff
  use error_data_module
  use discretization, only:nvar
  implicit none
  !
  type(crs_patcht),intent(in):: cg(0:dom_max_grids-1)
  type(sr_patcht),intent(in):: srg(-dom_max_sr_grids:0)
  type(error_data),intent(out)::errors(0:ncgrids-1+nvcycles+nsr)
  external Apply1, Apply2, Relax
  double precision,intent(out)::Res0 ! |f| on finest grid
  integer,intent(out):: nViters
  !     Local Vars
  integer:: iter,lev
  double precision:: Res,Reslast
  double precision,target,dimension(& ! alloc fine grid data here
       cg(0)%p%all%lo%i:cg(0)%p%all%hi%i,&
       cg(0)%p%all%lo%j:cg(0)%p%all%hi%j,&
       cg(0)%p%all%lo%k:cg(0)%p%all%hi%k,nvar)::L2u_f,f0,ux
  type(data_ptr),dimension(-nsr:0)::sr_ux,sr_f,sr_aux,sr_rhs
  interface
     subroutine MGSR(srg,cg,fl,sr_ux,sr_f,sr_aux,sr_rhs,uxC0,auxC0,Apply1,Apply2,Relax,Res0,errors)
       use base_data_module
       use pms, only:nvcycles,verbose,nsr,ncgrids,dom_max_sr_grids,dom_max_grids
       use error_data_module
       use mpistuff,only:mype
       use discretization, only:nvar
       implicit none
       type(sr_patcht),intent(in)::srg(-dom_max_sr_grids:0)
       type(crs_patcht),intent(in)::cg(0:dom_max_grids-1)
       type(data_ptr),dimension(-nsr:0)::sr_ux,sr_f,sr_aux,sr_rhs 
       double precision,intent(in),dimension(&
            cg(0)%p%all%lo%i:cg(0)%p%all%hi%i,& 
            cg(0)%p%all%lo%j:cg(0)%p%all%hi%j,& 
            cg(0)%p%all%lo%k:cg(0)%p%all%hi%k,nvar)::auxC0,uxC0
       type(error_data),intent(out)::errors(0:ncgrids-1+nvcycles+nsr)
       integer,intent(in)::fl
       double precision,intent(out)::Res0
       external Apply1, Apply2, Relax
     end subroutine MGSR
  end interface

  call formF(f0,cg(0)%p,cg(0)%p%max,cg(0)%t%ipe)
  Res0 = norm(f0,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
  if (verbose.gt.1.and.nsr.eq.0.and.nvcycles>0) then
     call formExactU(L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%t%ipe,cg(0)%p%dx)
     Reslast = norm(L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2) 
     if (mype==0) write(6,'(A,E12.4,A,E12.4,A,I5)') '0) solve: |f|_2=',&
          Res0,', |u_exact|_2=',Reslast,', nx=',(cg(0)%p%max%hi%i-cg(0)%p%max%lo%i+1)*cg(0)%t%npe%i
  end if
  Reslast = Res0
  ! FMG (SR coarse grid) solve
  iter = 0
  err_lev = 0
  if (nfcycles .ne. 0) then
     iter = iter + 1
     ! SR start of FMG
     call MGF(ux,f0,cg,0,Apply1,Apply2,Relax,errors)
     ! diagnostics
     if (verbose.gt.1.and.nsr==0.and.nvcycles>0) then
        call Apply2(L2u_f,ux,cg(0)%p,cg(0)%t)
        Res = norm(f0-L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
        if (mype==0)write(6,'(I2,A,E12.4,A,F10.6)') &
             iter,') solve: FMG done |f-Au|_2=',Res,', rate=',Res/Reslast
        Reslast=Res
     else if (nsr==0.and.nvcycles>0) then
        Res = norm(f0-L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
     end if
     ! SR 
     if (nsr.gt.0) then
        NULLIFY(sr_ux(0)%ptr,sr_f(0)%ptr,sr_aux(0)%ptr,sr_rhs(0)%ptr)
        do lev=-1,-nsr,-1
           NULLIFY(sr_ux(lev)%ptr,sr_f(lev)%ptr,sr_aux(lev)%ptr,sr_rhs(lev)%ptr)
           call MGSR(srg,cg,lev,sr_ux,sr_f,sr_aux,sr_rhs,ux,L2u_f,Apply1,Apply2,Relax,Res0,errors)
        end do
        
     end if
  else
     ! initailize what is done in FMG
     Res = Res0
     ux=0.d0
     if (nsr.ne.0) stop 'SR but no F-cycle'
  end if
 
  ! finish with V-cycle - no V with SR
  nViters = 0
  do while(nsr==0.and.iter<nvcycles.and.Res/Res0>rtol)
     call MGV(ux,f0,cg,0,Apply1,Relax)
     iter=iter+1 ! 2 for FMG and 1 for V 
     ! residual
     call Apply2(L2u_f,ux,cg(0)%p,cg(0)%t)            
     Res = norm(f0-L2u_f,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
     if (mype==0 .and. verbose.gt.0)write(6,'(I2,A,E12.4,A,F10.6)') &
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
        DEALLOCATE(sr_ux(lev)%ptr,sr_f(lev)%ptr,sr_aux(lev)%ptr,sr_rhs(lev)%ptr)
     end do
  end if
  return
end subroutine solve
!-----------------------------------------------------------------------
! MGF: recursive F-cycle with data on stack
!-----------------------------------------------------------------------
recursive subroutine MGF(ux,f,g,lev,Apply1,Apply2,Relax,errors)
  use grid_module
  use mpistuff
  use error_data_module
  use pms
  use discretization, only:nvar
  implicit none
  integer,intent(in) :: lev
  type(crs_patcht),intent(in):: g(0:dom_max_grids-1)
  type(error_data),intent(out)::errors(0:ncgrids-1+nvcycles+nsr)
  double precision,intent(in),dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,& 
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,& 
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k, nvar) :: f ! no tau
  double precision,intent(out),dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,&
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,& 
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k, nvar) :: ux
  external Apply1, Apply2, Relax
  ! Local Vars
  integer :: ii,jj,kk,ist,jst,kst
  double precision::tt
  double precision,dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,&
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,&
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k,nvar)::fTmp
  ! Coarse data. This is new data on the stack
  double precision,dimension(&
       g(lev+1)%p%all%lo%i:g(lev+1)%p%all%hi%i,&
       g(lev+1)%p%all%lo%j:g(lev+1)%p%all%hi%j,&
       g(lev+1)%p%all%lo%k:g(lev+1)%p%all%hi%k, nvar):: uxC,fC
  type(ipoint)::cfoffset

  cfoffset%i = 0;   cfoffset%j = 0;   cfoffset%k = 0; ! non-SR offsets

  ux=0.d0 ! zero out new level. everthing uses initial solution except FMG
  bc_valid = g(lev)%p%max%hi ! debug
  if ( is_top(g(lev)) ) then
     ! the real start of FMG - coarse grid solve
     if (lev.gt.0 .and. is_top(g(lev-1))) stop 'MGF: not coarsest grid'
     !if (mype==0.and.verbose.gt.2)write(6,'(A,I2)') '[ 0]MGF: bot solve lev=',nsr+lev
     !call Relax(ux,f,g(lev)%p,g(lev)%t,ncoarsesolveits)  ! coar grid solve
     call MGV(ux,f,g,lev,Apply1,Relax)
     if (verbose.gt.2) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'lev=',nsr+lev,') MGF: after bot solve |u|=',tt
     end if
  else
     ! go "down", allocating data (on stack), forming F, zero out u
     ! form F explicity on fine grid
     call formF(fC,g(lev+1)%p,g(lev+1)%p%max,g(lev+1)%t%ipe) 
     call MGF(uxC,fC,g,lev+1,Apply1,Apply2,Relax,errors) 
     if (verbose.gt.2) then
        tt = norm(uxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'lev=',nsr+lev+1,') pre FMG prol |u_c|=',tt
     end if
     ! FMG prolongate
     call Prolong(ux,uxC,g(lev)%t,g(lev+1)%t,g(lev)%p,g(lev+1)%p,cfoffset,.true.)
     
     if (verbose.gt.2) then
        tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'lev=',nsr+lev,  ') FMG prol dest |u|=',tt
     end if

     if (nsmoothsfmg>0) then
        ! pre FMG smoothing
        call Relax(ux,f,g(lev)%p,g(lev)%t,nsmoothsfmg)
        if (verbose.gt.2) then
           tt = norm(ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
           if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: after FMG pre sm |u|=',tt
        end if
     end if

     ! start of multiple v-cycle
     jj = nfmgvcycles !; if (lev==0) jj=10 ! accurate coarse solve
     do ii=1,jj
        call MGV(ux,f,g,lev,Apply1,Relax)
     enddo
  end if

  ! form errors & convergance measure, recursive so goes from coarse to fine (good)
  call formExactU(fTmp,g(lev)%p%all,g(lev)%p%max,g(lev)%t%ipe,g(lev)%p%dx)
  errors(err_lev)%uerror = norm(fTmp-ux,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,error_norm)
  call Apply2(fTmp,ux,g(lev)%p,g(lev)%t)
  tt=norm(fTmp-f,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
  errors(err_lev)%resid=tt
  ! like SR
  if (verbose.gt.2) then
     if (mype==0)write(6,'(A,I2,A,E12.4,A,E12.4,A,I8,I8,I8)') &
          '     lev=',nsr+lev,') MGF |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror,', n=',&
          g(lev)%p%max%hi%i*g(lev)%t%npe%i,&
          g(lev)%p%max%hi%j*g(lev)%t%npe%j,&
          g(lev)%p%max%hi%k*g(lev)%t%npe%k
  else if (verbose.gt.0) then
     if (mype==0)write(6,'(A,I2,A,E12.4,A,E12.4)') &
          '     lev=',nsr+lev,') MGF |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror
  end if
  err_lev = err_lev + 1

  return
end subroutine MGF
!-----------------------------------------------------------------------
recursive subroutine MGV(uxF,fF,g,lev,Apply,Relax)
  use grid_module
  use mpistuff
  use pms
  use discretization
  implicit none
  integer,intent(in)::lev
  type(crs_patcht),intent(in)::g(0:dom_max_grids-1)
  double precision,dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i, g(lev)%p%all%lo%j:g(lev)%p%all%hi%j, &
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k,nvar) :: uxF
  double precision,intent(in),dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i, g(lev)%p%all%lo%j:g(lev)%p%all%hi%j, &
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k,nvar) :: fF ! has tau when called recursively
  external Apply,Relax
  !
  type(ipoint)::cfoffset
  double precision::tt
  integer:: ii,jj,kk
  !     Coarse, here is the allocation on stack
  double precision,dimension(&
       g(lev+1)%p%all%lo%i:g(lev+1)%p%all%hi%i,&
       g(lev+1)%p%all%lo%j:g(lev+1)%p%all%hi%j, &
       g(lev+1)%p%all%lo%k:g(lev+1)%p%all%hi%k,nvar) :: uxC, tmpC, auxC
  ! This level
  double precision,dimension(&
       g(lev)%p%all%lo%i:g(lev)%p%all%hi%i,&
       g(lev)%p%all%lo%j:g(lev)%p%all%hi%j,& 
       g(lev)%p%all%lo%k:g(lev)%p%all%hi%k, nvar) :: auxF
  
  ! non-SR need a dummy for the offset
  cfoffset%i = 0; cfoffset%j = 0;   cfoffset%k = 0; 

  if ( is_top(g(lev)) ) then
     call Relax(uxf,fF,g(lev)%p,g(lev)%t,ncoarsesolveits)
     if (verbose.gt.2) then
        tt = norm(uxF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: after bot solve |u|=',tt
     end if
  else
     if (verbose.gt.2) then
        tt = norm(fF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')   '       lev=',nsr+lev,') V: |f_0|=',tt
        tt = norm(uxF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4,I4)')'       lev=',nsr+lev,') V: |u_0|=',tt,g(lev)%p%max%hi%i
     end if
     !     pre smoothing
     call Relax(uxF,fF,g(lev)%p,g(lev)%t,nsmoothsdown)
     if (verbose.gt.2) then
        tt = norm(uxF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: after pre sm |u|=',tt
     end if
     !     resid: f - Au -- f has tau in it
     call Apply(auxF,uxF,g(lev)%p,g(lev)%t)
     auxF = fF - auxF

     if (verbose.gt.2) then
        tt = norm(auxF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: |f-Au|=',tt
        tt = norm(fF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: |f|=',tt
     end if
     bc_valid = g(lev+1)%p%max%hi ! debug
     !     restrict solution & f-Au
     call RestrictFuse(g(lev+1)%p,g(lev)%p,g(lev+1)%t,cfoffset,uxC,auxC,uxF,auxF)
     if (verbose.gt.2) then
        tt = norm(auxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev+1,') V: after R |f-Au|=',tt
        tt = norm(uxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev+1,') V: after R |u|=',tt
     end if
     !     rhsC = Fc + A*R(u) - R(A*u) = Fc + tau_c
     call Apply(tmpC,uxC,g(lev+1)%p,g(lev+1)%t)
     !call formF(fC,g(lev+1)%p,g(lev+1)%p%max,g(lev+1)%t%ipe) ! this is wasteful, could cache
     !rhsC = fC + tmpC - rhsC
     auxC = auxC + tmpC
     if (verbose.gt.2) then
        tt = norm(auxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev+1,') V: |f+tau_c|=',tt
        tt = norm(tmpC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev+1,') V: |AU_c|=',tt
     end if
     !     Temporarily store uxC into tmpC
     tmpC = uxC
     !     start of v(1) or w(2) cycle 
     jj = ncycles ! ; if (lev==0) jj=40 ! two level method, sort of
     do ii=1,jj
        call MGV(uxC,auxC,g,lev+1,Apply,Relax)
     enddo
     if (verbose.gt.2) then
        tt = norm(uxC,g(lev+1)%p%all,g(lev+1)%p%max,g(lev+1)%p%dx,g(lev+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev+1,') V: after MGV() |u|=',tt
     end if
     !     subtract off old Uc
     tmpC = uxC - tmpC
     ! prolongate and correct
bc_valid = g(lev)%p%max%hi ! debug
     call Prolong(uxF,tmpC,g(lev)%t,g(lev+1)%t,g(lev)%p,g(lev+1)%p,cfoffset,.false.)
     if (verbose.gt.2) then
        tt = norm(uxF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: after prol |u|=',tt
     end if
     ! end of v cycle, post smoothing
     call Relax(uxF,fF,g(lev)%p,g(lev)%t,nsmoothsup)
     if (verbose.gt.2) then
        tt = norm(uxF,g(lev)%p%all,g(lev)%p%max,g(lev)%p%dx,g(lev)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: after post |u|=',tt
     end if
  end if
  return
end subroutine MGV
!-----------------------------------------------------------------------
! SR: copy_sr_crs2 
!  - copy two vectors (u,f) from SR (sr0) valid to normal (cg0) 'max'
!-----------------------------------------------------------------------
subroutine copy_sr_crs2(crs1,crs2,sr1,sr2,cg0,srg0)
  use bc_module
  use pms
  use mpistuff,only:mype
  use discretization,only:nvar,nsg
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
  !
  integer::ii,jj,kk,sri,srj,srk
  ! copy
  ii=0
  do sri=srg0%val%lo%i,srg0%val%hi%i
     ii=ii+1
     jj=0
     do srj=srg0%val%lo%j,srg0%val%hi%j
        jj=jj+1
        kk=0
        do srk=srg0%val%lo%k,srg0%val%hi%k
           kk=kk+1
           crs1(ii,jj,kk,:) = sr1(sri,srj,srk,:)
           crs2(ii,jj,kk,:) = sr2(sri,srj,srk,:)
        end do
     end do
  end do
  call SetBCs(crs1,cg0%p,cg0%t,nsg) ! normal BC exchange 
  return
end subroutine copy_sr_crs2
!
subroutine sr_exchange(u_sr,srg0,cg0)
  use base_data_module
  use mpistuff,only:mype
  use discretization
  implicit none
  type(sr_patcht),intent(in)::srg0
  type(crs_patcht),intent(in)::cg0
  double precision,dimension(&
       srg0%p%all%lo%i:srg0%p%all%hi%i,&
       srg0%p%all%lo%j:srg0%p%all%hi%j,&
       srg0%p%all%lo%k:srg0%p%all%hi%k,nvar)::u_sr
  !
  type(patcht)::p
  type(ipoint)::ng
  type(ipoint)::shift
  ng = nsg
  shift%i=srg0%val%lo%i-1
  shift%j=srg0%val%lo%j-1
  shift%k=srg0%val%lo%k-1
  ! set ghost sizes
  if (srg0%p%max%hi%i-srg0%val%hi%i+nsg%i .gt. ng%i) ng%i = srg0%p%max%hi%i-srg0%val%hi%i+nsg%i
  if (srg0%val%lo%i-1+nsg%i .gt. ng%i) ng%i = srg0%val%lo%i-1+nsg%i
  ! 
  if (srg0%p%max%hi%j-srg0%val%hi%j+nsg%j .gt. ng%j) ng%j = srg0%p%max%hi%j-srg0%val%hi%j+nsg%j
  if (srg0%val%lo%j-1+nsg%j .gt. ng%j) ng%j = srg0%val%lo%j-1+nsg%j
  ! 
  if (srg0%p%max%hi%k-srg0%val%hi%k+nsg%k .gt. ng%k) ng%k = srg0%p%max%hi%k-srg0%val%hi%k+nsg%k
  if (srg0%val%lo%k-1+nsg%k .gt. ng%k) ng%k = srg0%val%lo%k-1+nsg%k
  if (nsg%i.gt.srg0%p%max%hi%i+nsg%i) then
     print *,ng%i,srg0%val%hi%i-srg0%val%lo%i+1,cg0%p%all%hi%i-cg0%p%all%lo%i+1
     stop 'copy_crs_sr_w_bc_ex SR buffer too big for (valid) domain. increase box size'
  end if
  ! strange patch with big exchange areas
  p%dx%i = 0.d0 ! not used
  p%dx%j = 0.d0 
  p%dx%k = 0.d0 
  ! shift alloc region to put valid in max
  p%all%lo%i=srg0%p%all%lo%i-shift%i
  p%all%hi%i=srg0%p%all%hi%i-shift%i
  p%all%lo%j=srg0%p%all%lo%j-shift%j
  p%all%hi%j=srg0%p%all%hi%j-shift%j
#ifndef TWO_D
  p%all%lo%k=srg0%p%all%lo%k-shift%k
  p%all%hi%k=srg0%p%all%hi%k-shift%k
#else
  p%all%lo%k=1
  p%all%hi%k=1
#endif 
  ! max is valid size
  p%max%hi%i=srg0%val%hi%i-srg0%val%lo%i+1
  p%max%hi%j=srg0%val%hi%j-srg0%val%lo%j+1
#ifndef TWO_D
  p%max%hi%k=srg0%val%hi%k-srg0%val%lo%k+1
#else
  p%max%hi%k=1
#endif
  call SetBCs(u_sr,p,cg0%t,ng)

  return
end subroutine sr_exchange
!-----------------------------------------------------------------------
! SR: copy_crs_sr_w_bc_ex - copy vector from normal to SR
!-----------------------------------------------------------------------
subroutine copy_crs_sr_w_bc_ex(u_sr,u_crs,srg0,cg0)
  use mpistuff,only:mype
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
  !
  integer::ii,jj,kk,sri,srj,srk,nbuff
  ! copy to valid from whole
  ii=0
  do sri=srg0%val%lo%i,srg0%val%hi%i
     ii=ii+1
     jj=0
     do srj=srg0%val%lo%j,srg0%val%hi%j
        jj=jj+1
        kk=0
        do srk=srg0%val%lo%k,srg0%val%hi%k
           kk=kk+1
           u_sr(sri,srj,srk,:) = u_crs(ii,jj,kk,:)
        end do
     end do
  end do
  if (ii/=cg0%p%max%hi%i) stop 'copy_crs_sr_w_bc_ex ii'
  if (jj/=cg0%p%max%hi%j) stop 'copy_crs_sr_w_bc_ex jj'
  if (kk/=cg0%p%max%hi%k) stop 'copy_crs_sr_w_bc_ex kk'

  call sr_exchange(u_sr,srg0,cg0)

  return
end subroutine copy_crs_sr_w_bc_ex
!-----------------------------------------------------------------------
! SR: MSRG - SR part of FMG
!-----------------------------------------------------------------------
subroutine MGSR(srg,cg,lev,sr_ux,sr_f,sr_aux,sr_rhs,uxC0,auxC0,Apply1,Apply2,Relax,Res0,errors)
  use grid_module
  use pms
  use error_data_module
  use mpistuff,only:mype
  use sr_interfaces
  implicit none
  type(sr_patcht),intent(in)::srg(-dom_max_sr_grids:0)
  type(crs_patcht),intent(in)::cg(0:dom_max_grids-1)
  type(data_ptr),dimension(-nsr:0)::sr_ux,sr_f,sr_aux,sr_rhs
  double precision,dimension(&
       cg(0)%p%all%lo%i:cg(0)%p%all%hi%i,& 
       cg(0)%p%all%lo%j:cg(0)%p%all%hi%j,& 
       cg(0)%p%all%lo%k:cg(0)%p%all%hi%k,nvar)::uxC0,auxC0 ! in/out, buffer
  type(error_data),intent(out)::errors(0:ncgrids-1+nvcycles+nsr)
  integer,intent(in)::lev ! negative index into srg(lev)
  double precision,intent(out)::Res0 ! residual norm on finest grid
  external Apply1, Apply2, Relax
  ! 
  double precision,pointer,DIMENSION(:,:,:,:)::uxF,auxF,rhsF,uxC,rhsC,auxC,tmpC,fTmp
  type(data_ptr),dimension(lev+1:0)::tmpSRC
  integer::iflv,AllocateStatus,ii,jj,kk
  double precision::tt,t2

  ! allocate coarse grid if needed, first call
  if (.not.ASSOCIATED(sr_ux(0)%ptr)) then
     if (lev.ne.-1) stop 'MGSR: lev.ne.-1????'
     ALLOCATE(sr_ux(0)%ptr(&
       srg(0)%p%all%lo%i:srg(0)%p%all%hi%i,&
       srg(0)%p%all%lo%j:srg(0)%p%all%hi%j,&
       srg(0)%p%all%lo%k:srg(0)%p%all%hi%k,nvar),sr_f(0)%ptr(&
       srg(0)%p%all%lo%i:srg(0)%p%all%hi%i,&
       srg(0)%p%all%lo%j:srg(0)%p%all%hi%j,&
       srg(0)%p%all%lo%k:srg(0)%p%all%hi%k,nvar),sr_aux(0)%ptr(&
       srg(0)%p%all%lo%i:srg(0)%p%all%hi%i,&
       srg(0)%p%all%lo%j:srg(0)%p%all%hi%j,&
       srg(0)%p%all%lo%k:srg(0)%p%all%hi%k,nvar),sr_rhs(0)%ptr(&
       srg(0)%p%all%lo%i:srg(0)%p%all%hi%i,&
       srg(0)%p%all%lo%j:srg(0)%p%all%hi%j,&
       srg(0)%p%all%lo%k:srg(0)%p%all%hi%k,nvar), STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "MGSR: *** Not enough memory ***"
     ! copy in solution, apply BC, and exchange
     uxC => sr_ux(0)%ptr 
bc_valid = srg(0)%val%hi ! debug
     call copy_crs_sr_w_bc_ex(uxC,uxC0,srg(0),cg(0)) ! input
     
     tt = norm(uxC0,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
!!$     if (mype==0)write(6,'(A,I2,A,E12.4)')        'lev=',nsr,') MGSR copy_crs_sr_w_bc_ex |u_src|=',tt
     t2 = norm(uxC,srg(0)%p%all,srg(0)%val,srg(0)%p%dx,srg(0)%t%comm,2)
!!$     if (mype==0)write(6,'(A,I2,A,E12.4,A,E12.4)')'lev=',nsr,') MSRG copy_crs_sr_w_bc_ex |u_dest|=',t2,', error:',tt-t2
     if (abs(tt-t2)/t2>1.d-12) stop 'MSRG: big error'



     ! just compute F, could copy from normal grid 0
     fTmp => sr_f(0)%ptr
     call formF(fTmp,srg(0)%p,srg(0)%val,srg(0)%t%ipe)
  end if

  ! allocate new fine grid
  ALLOCATE(sr_ux(lev)%ptr(&
       srg(lev)%p%all%lo%i:srg(lev)%p%all%hi%i,&
       srg(lev)%p%all%lo%j:srg(lev)%p%all%hi%j,&
       srg(lev)%p%all%lo%k:srg(lev)%p%all%hi%k,nvar),sr_f(lev)%ptr(&
       srg(lev)%p%all%lo%i:srg(lev)%p%all%hi%i,&
       srg(lev)%p%all%lo%j:srg(lev)%p%all%hi%j,&
       srg(lev)%p%all%lo%k:srg(lev)%p%all%hi%k,nvar),sr_aux(lev)%ptr(&
       srg(lev)%p%all%lo%i:srg(lev)%p%all%hi%i,&
       srg(lev)%p%all%lo%j:srg(lev)%p%all%hi%j,&
       srg(lev)%p%all%lo%k:srg(lev)%p%all%hi%k,nvar),sr_rhs(lev)%ptr(&
       srg(lev)%p%all%lo%i:srg(lev)%p%all%hi%i,&
       srg(lev)%p%all%lo%j:srg(lev)%p%all%hi%j,&
       srg(lev)%p%all%lo%k:srg(lev)%p%all%hi%k,nvar),STAT=AllocateStatus)
  if (allocatestatus /= 0) stop "MGSR: *** not enough memory ***"

  ! SR1 = Prol + SR2
  fTmp => sr_f(lev)%ptr ! form F on new level and use for first relax
  call formF(fTmp,srg(lev)%p,srg(lev)%val,srg(lev)%t%ipe)
  Res0 = norm(fTmp,srg(lev)%p%all,srg(lev)%val,srg(lev)%p%dx,srg(lev)%t%comm,2) ! out
  if (verbose.gt.2) then
     if (mype==0) print *,'SR:'
     uxC => sr_ux(lev+1)%ptr 
     tt = norm(uxC,srg(lev+1)%p%all,srg(lev+1)%val,srg(lev+1)%p%dx,srg(lev+1)%t%comm,2) 
     if (mype==0)write(6,'(A,I2,A,E12.4)')'lev=',nsr+lev+1,') pre FMG prol |u_c|=',tt
  end if
  bc_valid = srg(lev)%val%hi ! debug

  ! FMG prolongation
  uxF => sr_ux(lev)%ptr 
  uxF = 0.d0
  uxC => sr_ux(lev+1)%ptr 
  call Prolong(uxF,uxC,srg(lev)%t,srg(lev+1)%t,srg(lev)%p,srg(lev+1)%p,srg(lev+1)%cfoffset,.true.)
  if (verbose.gt.2) then
     if (mype==0) print *,'SR:'
     tt = norm(uxF,srg(lev)%p%all,srg(lev)%val,srg(lev)%p%dx,srg(lev)%t%comm,2)
     if (mype==0)write(6,'(A,I2,A,E12.4)')'lev=',nsr+lev,  ') FMG prol dest |u|=',tt
  end if
  
  if (nsmoothsfmg>0) then
     uxF  => sr_ux (lev)%ptr
     fTmp => sr_f(lev)%ptr ! no SR
     ! SR2:    pre smoothing
     call Relax(uxF,fTmp,srg(lev)%p,srg(lev)%t,nsmoothsfmg)
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(uxF,srg(lev)%p%all,srg(lev)%val,srg(lev)%p%dx,srg(lev)%t%comm,2)
        if (mype==0) write(6,'(A,I2,A,E12.4)')'       lev=',nsr+lev,') V: after FMG pre sm |u|=',tt
     end if
  end if

  ! SR2 = smooth restrict u & f-Au down to coarsest grid 0 (lev+1)
  fTmp => sr_f(lev)%ptr
  rhsF => sr_rhs(lev)%ptr
  rhsF = fTmp ! first (finest) level has no tau correction, just F
  do iflv=lev,-1,1
     uxF  => sr_ux (iflv)%ptr
     !fF => sr_f(iflv)%ptr
     auxF => sr_aux(iflv)%ptr
     rhsF => sr_rhs(iflv)%ptr
     uxC  => sr_ux (iflv+1)%ptr ! level 0 in last loop
     !fC => sr_f(iflv+1)%ptr
     auxC => sr_aux(iflv+1)%ptr
     rhsC => sr_rhs(iflv+1)%ptr
     ALLOCATE(tmpSRC(iflv+1)%ptr(&
          srg(iflv+1)%p%all%lo%i:srg(iflv+1)%p%all%hi%i,&
          srg(iflv+1)%p%all%lo%j:srg(iflv+1)%p%all%hi%j,&
          srg(iflv+1)%p%all%lo%k:srg(iflv+1)%p%all%hi%k,nvar),STAT=AllocateStatus)
     IF (AllocateStatus /= 0) STOP "MSRG: *** Not enough memory ***"
     tmpC => tmpSRC(iflv+1)%ptr
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(rhsF,srg(iflv)%p%all,srg(iflv)%val,srg(iflv)%p%dx,srg(iflv)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')   '       lev=',nsr+iflv,') V: |f_0|=',tt
        tt = norm(uxF,srg(iflv)%p%all,srg(iflv)%val,srg(iflv)%p%dx,srg(iflv)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4,I4)')'       lev=',nsr+iflv,') V: |u_0|=',tt,srg(iflv)%val%hi%i
     end if
     ! SR2:    pre smoothing
     call Relax(uxF,rhsF,srg(iflv)%p,srg(iflv)%t,nsmoothsdown)
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(uxF,srg(iflv)%p%all,srg(iflv)%val,srg(iflv)%p%dx,srg(iflv)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv,') V: after pre sm |u|=',tt
     end if
     !     form f-Au ! f has tau ?? 
     call Apply1(auxF,uxF,srg(iflv)%p,srg(iflv)%t)
     auxF = rhsF - auxF
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(auxF,srg(iflv)%p%all,srg(iflv)%val,srg(iflv)%p%dx,srg(iflv)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv,') V: |f-Au|=',tt
        tt = norm(rhsF,srg(iflv)%p%all,srg(iflv)%val,srg(iflv)%p%dx,srg(iflv)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv,') V: |f|=',tt
     end if
     bc_valid = srg(iflv+1)%val%hi ! debug
     !     restrict u & f-Au
     !if (nsmoothsfmg.le.0.and.nsmoothsdown.le.0.and.iflv==lev) tmpC = uxC ! lets keep the solution if no smoothing
     call RestrictFuse(srg(iflv+1)%p,srg(iflv)%p,srg(iflv+1)%t,srg(iflv+1)%cfoffset,uxC,auxC,uxF,auxF)
     !if (nsmoothsfmg.le.0.and.nsmoothsdown.le.0.and.iflv==lev) uxC = tmpC ! lets keep the solution if no smoothing
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(auxC,srg(iflv+1)%p%all,srg(iflv+1)%val,srg(iflv+1)%p%dx,srg(iflv+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv+1,') V: after R |f-Au|=',tt
        tt = norm(uxC,srg(iflv+1)%p%all,srg(iflv+1)%val,srg(iflv+1)%p%dx,srg(iflv+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv+1,') V: after R |u|=',tt
     end if
     !     rhsC = Fc + A*R(u) - R(A*u) = Fc + tau_c
     call Apply1(tmpC,uxC,srg(iflv+1)%p,srg(iflv+1)%t)
     rhsC = auxC + tmpC  
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(rhsC,srg(iflv+1)%p%all,srg(iflv+1)%val,srg(iflv+1)%p%dx,srg(iflv+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv+1,') V: |f+tau_c|=',tt
        tt = norm(tmpC,srg(iflv+1)%p%all,srg(iflv+1)%val,srg(iflv+1)%p%dx,&
             srg(iflv+1)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv+1,') V: |AU_c|=',tt
     end if
     !     Temporarily copy uxC into tmpSRC
     tmpC = uxC 
  end do

  ! get U & F from (these are already set)
  rhsC => sr_rhs(0)%ptr ! rhsC in MGV
  uxC  => sr_ux (0)%ptr 
  ! 'copy out' to coarse grid
  call copy_sr_crs2(uxC0,auxC0,uxC,rhsC,cg(0),srg(0))

  ! debug
  !tt = norm(auxC0,cg(0)%p%all,cg(0)%p%max,cg(0)%p%dx,cg(0)%t%comm,2)
  !if (mype==0)write(6,'(A,I2,A,E12.4)')'    lev=',0,') SR: |f_crs_0|=',tt
  !t2 = norm(fC,srg(0)%p%all,srg(0)%val,srg(0)%p%dx,srg(0)%t%comm,2)
  !if (mype==0)write(6,'(A,I2,A,E12.4)')'    lev=',0,') SR: |f_sr_0|=',t2
  !if (abs(tt-t2)/t2>1.d-12) stop 'MSRG: big error (2)'

  ! coarse grid solve - these pointers are coarse data types, no SR buffers
  call MGV(uxC0,auxC0,cg,0,Apply1,Relax)

  ! 'copy in' coarse grid do SR coarse grid (cg.0-->sr.0)  
  call copy_crs_sr_w_bc_ex(uxC,uxC0,srg(0),cg(0))
  if (verbose.gt.2) then
     if (mype==0) print *,'SR:'
     tt = norm(uxC,srg(0)%p%all,srg(0)%val,srg(0)%p%dx,srg(0)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr,') V: after MGV() |u|=',tt
  end if

  ! SR3: prolongate, update, smooth, coarse to fine
  do iflv=-1,lev,-1
     uxF  =>  sr_ux(iflv)%ptr
     !fF => sr_f(iflv)%ptr
     auxF => sr_aux(iflv)%ptr
     rhsF => sr_rhs(iflv)%ptr
     uxC => sr_ux(iflv+1)%ptr
     !fC => sr_f(iflv+1)%ptr
     auxC => sr_aux(iflv+1)%ptr
     rhsC => sr_rhs(iflv+1)%ptr
     tmpC => tmpSRC(iflv+1)%ptr
     !     subtract off old Uc
     tmpC = uxC - tmpC ! make increment now
     ! prolongate and correct
     bc_valid = srg(iflv)%val%hi ! debug
     ! u = u + P(u_c - u_c_0)
     call Prolong(uxF,tmpC,srg(iflv)%t,srg(iflv+1)%t,&
          srg(iflv)%p,srg(iflv+1)%p,srg(iflv+1)%cfoffset,.false.)
     DEALLOCATE(tmpC)
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(uxF,srg(iflv)%p%all,srg(iflv)%val,srg(iflv)%p%dx,srg(iflv)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv,') V: after prol |u|=',tt
     end if
     ! end of v cycle, post smoothing
     call Relax(uxF,rhsF,srg(iflv)%p,srg(iflv)%t,nsmoothsup)
     if (verbose.gt.2) then
        if (mype==0) print *,'SR:'
        tt = norm(uxF,srg(iflv)%p%all,srg(iflv)%val,srg(iflv)%p%dx,srg(iflv)%t%comm,2)
        if (mype==0)write(6,'(A,I2,A,E12.4)')'       lev=',nsr+iflv,') V: after post |u|=',tt
     end if
  end do

  uxF  =>  sr_ux(lev)%ptr
  auxF => sr_aux(lev)%ptr
  fTmp => sr_f(lev)%ptr
  ! form errors & convergance measure, recursive so goes from coarse to fine (good)

  call formExactU(auxF,srg(lev)%p%all,srg(lev)%val,srg(lev)%t%ipe,srg(lev)%p%dx)
  tt=norm(auxF-uxF,srg(lev)%p%all,srg(lev)%val,srg(lev)%p%dx,srg(lev)%t%comm,error_norm)
  errors(err_lev)%uerror=tt
  call sr_exchange(uxF,srg(lev),cg(0)) ! get an accurate residual
  call Apply2(auxF,uxF,srg(lev)%p,srg(lev)%t) 
  tt=norm(fTmp-auxF,srg(lev)%p%all,srg(lev)%val,srg(lev)%p%dx,srg(lev)%t%comm,2)
  errors(err_lev)%resid=tt
  ! like FMG
  if (verbose.gt.2) then
     if (mype==0) print *,'SR:'
     if (mype==0)write(6,'(A,I2,A,E12.4,A,E12.4,A,I8,I8,I8)') &
          '     lev=',nsr+lev,') MGF |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror,', n=',&
          (srg(lev)%val%hi%i-srg(lev)%val%lo%i+1)*srg(lev)%t%npe%i,&
          (srg(lev)%val%hi%j-srg(lev)%val%lo%j+1)*srg(lev)%t%npe%j,&
          (srg(lev)%val%hi%k-srg(lev)%val%lo%k+1)*srg(lev)%t%npe%k
  else if (verbose.gt.0) then
     if (mype==0)write(6,'(A,I2,A,E12.4,A,E12.4)') &
          '     lev=',nsr+lev,') MGF-SR |res|_2=',errors(err_lev)%resid,&
          ', |error|=',errors(err_lev)%uerror
  end if
  err_lev = err_lev + 1
  
  return
end subroutine MGSR
