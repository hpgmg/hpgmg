!-----------------------------------------------------------------
!     Ravi Samtaney & Mark Adams
!     Copyright 2014
!-----------------------------------------------------------------

!-----------------------------------------------------------------------
! RestrictFuse: fuse two restricts (u,f).
!-----------------------------------------------------------------------
subroutine RestrictFuse(cp,fp,ct,cOffset,uC,uC2,ux,ux2)
  use discretization,only:nvar,nsg
  use bc_module
  use mpistuff
  use tags
  use pms, only:mg_ref_ratio
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
  ! buffers are 1/2 of coarse
#ifdef TWO_D
  double precision,dimension(fp%max%hi%i/2,fp%max%hi%j/2,1,2*nvar)::&
       sb,rb1,rb2,rb3,rb0
#else
  double precision,dimension(fp%max%hi%i/2,fp%max%hi%j/2,fp%max%hi%k/2,2*nvar)::&
       sb,rb0,rb1,rb2,rb3,rb4,rb5,rb6,rb7
#endif
  double precision :: a11,a22,a33,a44,tmp
  integer::ic,jc,kc,pp,lnp
  integer::msg_id_recv(0:7),bufsz,msg_id_send(0:7)
  integer::ii,jj,kk
  type(ipoint)::csz,cstrt2,cend2,cend1

  if (mg_ref_ratio .ne. 2) stop 'RestrictFuse: not done for refratio != 2'

#ifdef HAVE_PETSC
  flops = 2*fp%max%hi%i*fp%max%hi%j*fp%max%hi%k
  call PetscLogEventBegin(events(3),ierr)
  call PetscLogFlops(flops,ierr)
#endif

  ! size of coarse patch to set & end of buffer
  csz%i = fp%max%hi%i/2
  csz%j = fp%max%hi%j/2
  csz%k = fp%max%hi%k/2
#ifdef TWO_D
  a11 = 9.d0/64.d0
  a22 = 3.d0/64.d0
  a33 = 1.d0/64.d0
#else
  a11 = 27.d0/64.d0/4.d0
  a22 = 9.d0/64.d0/4.d0
  a33 = 3.d0/64.d0/4.d0
  a44 = 1.d0/64.d0/4.d0 ! not used
#endif

  ! split grid
  if (fp%max%hi%i==cp%max%hi%i) then ! could this happen on an SR grid?
     if (fp%max%hi%j.ne.cp%max%hi%j .or. fp%max%hi%k.ne.cp%max%hi%k) stop 'RestrictFuse'
     ! start of 2nd half
     cstrt2%i = 1+csz%i
     cstrt2%j = 1+csz%j
     cstrt2%k = 1+csz%k
     ! end of 1st half
     cend1%i = csz%i
     cend1%j = csz%j
     cend1%k = csz%k
     ! end of 2nd half
     cend2%i = cstrt2%i+csz%i-1
     cend2%j = cstrt2%j+csz%j-1
     cend2%k = cstrt2%k+csz%k-1
#ifdef TWO_D
     bufsz = 2*nvar*csz%j*csz%i
     lnp = 4
#else 
     bufsz = 2*nvar*csz%j*csz%i*csz%k
     lnp = 8
#endif
     ! post recieves
     call MPI_Irecv(rb0,bufsz,MPI_DOUBLE_PRECISION,&
          0,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(0),ierr)
     call MPI_Irecv(rb1,bufsz,MPI_DOUBLE_PRECISION,&
          1,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(1),ierr)
     call MPI_Irecv(rb2,bufsz,MPI_DOUBLE_PRECISION,&
          2,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(2),ierr)
     call MPI_Irecv(rb3,bufsz,MPI_DOUBLE_PRECISION,&
          3,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(3),ierr)
#ifndef TWO_D  
     call MPI_Irecv(rb4,bufsz,MPI_DOUBLE_PRECISION,&
          4,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(4),ierr)
     call MPI_Irecv(rb5,bufsz,MPI_DOUBLE_PRECISION,&
          5,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(5),ierr)
     call MPI_Irecv(rb6,bufsz,MPI_DOUBLE_PRECISION,&
          6,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(6),ierr)
     call MPI_Irecv(rb7,bufsz,MPI_DOUBLE_PRECISION,&
          7,MSG_RESTRICT_TAG,ct%loc_comm,msg_id_recv(7),ierr)
#endif

     ! copy in restriction of fine patch into buffer
     kk=-1
     do kc=1,csz%k
        kk=kk+2 ! 1
        jj=-1
        do jc=1,csz%j
           jj=jj+2
           ii=-1
           do ic=1,csz%i
              ii=ii+2
#ifdef TWO_D
                 sb(ic,jc,kc,1:nvar)=.25d0*(&
                       ux(ii,  jj,  kk,:)+ux(ii+1,jj,  kk,:)&
                      +ux(ii+1,jj+1,kk,:)+ux(ii,  jj+1,kk,:))
                 sb(ic,jc,kc,nvar+1:2*nvar)=.25d0*(&
                      ux2(ii,   jj,  kk,:)+ux2(ii+1,jj,  kk,:)&
                      +ux2(ii+1,jj+1,kk,:)+ux2(ii,  jj+1,kk,:))
#else 
                 sb(ic,jc,kc,1:nvar)=.125d0*(&
                       ux(ii,  jj,  kk,  :)+ux(ii+1,jj,  kk,  :)&
                      +ux(ii+1,jj+1,kk,  :)+ux(ii,  jj+1,kk,  :)&
                      +ux(ii,  jj,  kk+1,:)+ux(ii+1,jj,  kk+1,:)&
                      +ux(ii+1,jj+1,kk+1,:)+ux(ii,  jj+1,kk+1,:))
                 sb(ic,jc,kc,nvar+1:2*nvar)=.125d0*(&
                       ux2(ii,  jj,  kk,  :)+ux2(ii+1,jj,  kk,  :)&
                      +ux2(ii+1,jj+1,kk,  :)+ux2(ii,  jj+1,kk,  :)&
                      +ux2(ii,  jj,  kk+1,:)+ux2(ii+1,jj,  kk+1,:)&
                      +ux2(ii+1,jj+1,kk+1,:)+ux2(ii,  jj+1,kk+1,:))
#endif
           enddo
        enddo
     enddo
     
     ! send to my cohort & me
     do pp=0,lnp-1            
        call MPI_Isend(sb,bufsz,MPI_DOUBLE_PRECISION,pp,&
             MSG_RESTRICT_TAG,ct%loc_comm,msg_id_send(pp),ierr)
     enddo
     MSG_RESTRICT_TAG = MSG_RESTRICT_TAG + 1

     !     0    recv (includes self recv)  (1,1[,1])
     call MPI_Wait(msg_id_recv(0), status, ierr)
     kk=1
     do kc=1,cend1%i
        jj=1
        do jc=1,cend1%j
           ii=1
           do ic=1,cend1%k
              uC (ic,jc,kc,:) = rb0(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb0(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
     !     1
     call MPI_Wait(msg_id_recv(1),status,ierr)
#ifdef TWO_D    
     jj=1         ! (1,2)
     do jc=cstrt2%j,cend2%j
        ii=1
        do ic=1,cend1%i
           uC (ic,jc,1,:) = rb1(ii,jj,1,1:nvar)
           uC2(ic,jc,1,:) = rb1(ii,jj,1,nvar+1:2*nvar)
           ii=ii+1
        enddo
        jj=jj+1
     enddo
#else 
     kk=1
     do kc=cstrt2%k,cend2%k   ! (1,1,2)
        jj=1
        do jc=1,cend1%j
           ii=1
           do ic=1,cend1%i
              uC (ic,jc,kc,:) = rb1(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb1(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
#endif
     !     2  
     call MPI_Wait(msg_id_recv(2), status, ierr)
#ifdef TWO_D
     jj=1
     do jc=1,cend1%j       ! (2,1)
        ii = 1
        do ic=cstrt2%i,cend2%i
           uC (ic,jc,kc,:) = rb2(ii,jj,1,1:nvar)
           uC2(ic,jc,kc,:) = rb2(ii,jj,1,nvar+1:2*nvar)
           ii = ii + 1
        enddo
        jj=jj+1
     enddo
#else 
     kk=1
     do kc=1,cend1%k       ! (1,2,1)
        jj=1 
        do jc=cstrt2%j,cend2%j
           ii=1
           do ic=1,cend1%i
               uC(ic,jc,kc,:) = rb2(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb2(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
#endif
     !     3
     call MPI_Wait(msg_id_recv(3), status, ierr)
#ifdef TWO_D
     jj=1                 ! (2,2)
     do jc=cstrt2%j,cend2%j
        ii=1
        do ic=cstrt2%i,cend2%i
           uC (ic,jc,1,:) = rb3(ii,jj,kk,1:nvar)
           uC2(ic,jc,1,:) = rb3(ii,jj,kk,nvar+1:2*nvar)
           ii=ii+1
        enddo
        jj=jj+1
     enddo
#else
     kk=1
     do kc=cstrt2%k,cend2%k  ! (1,2,2)
        jj=1
        do jc=cstrt2%j,cend2%j
           ii=1
           do ic=1,cend1%i
              uC (ic,jc,kc,:) = rb3(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb3(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
     !     4 (2,1,1)
     call MPI_Wait(msg_id_recv(4), status, ierr)
     kk=1
     do kc=1,cend1%k
        jj=1
        do jc=1,cend1%j
           ii=1
           do ic=cstrt2%i,cend2%i
              uC (ic,jc,kc,:) = rb4(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb4(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo     
        kk=kk+1
     enddo
     !     5 (2,1,2)
     call MPI_Wait(msg_id_recv(5), status, ierr)
     kk=1
     do kc=cstrt2%k,cend2%k 
        jj=1
        do jc=1,cend1%j
           ii=1
           do ic=cstrt2%i,cend2%i
              uC (ic,jc,kc,:) = rb5(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb5(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
     !     6 (2,2,1)
     call MPI_Wait(msg_id_recv(6), status, ierr)
     kk=1
     do kc=1,cend1%k
        jj=1
        do jc=cstrt2%j,cend2%j
           ii=1
           do ic=cstrt2%i,cend2%i
              uC (ic,jc,kc,:) = rb6(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb6(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
     !     7 (2,2,2)
     call MPI_Wait(msg_id_recv(7), status, ierr)
     kk=1
     do kc=cstrt2%k,cend2%k
        jj=1
        do jc=cstrt2%j,cend2%j
           ii=1
           do ic=cstrt2%i,cend2%i
              uC (ic,jc,kc,:) = rb7(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb7(ii,jj,kk,nvar+1:2*nvar)
              ii=ii+1
           enddo
           jj=jj+1
        enddo
        kk=kk+1
     enddo
#endif     
     do pp=0,lnp-1
        call MPI_Wait(msg_id_send(pp), status, ierr)
     enddo
  else ! normal grid (not split)
     ! c.max : |  offset  | | | | | offset |
!!$     write(6,'(I4,A,I4,I4,I4,A,I4,I4,I4,A,I4,I4,I4)') mype,&
!!$          '] fine max:',fp%max%hi%i,fp%max%hi%j,fp%max%hi%k,&
!!$          ', coarse offset: ',cOffset%i,cOffset%j,cOffset%k,&
!!$          ', coarse max:',cp%max%hi%i,cp%max%hi%j,cp%max%hi%k
     ! end of 1st half
     cend1%i = cOffset%i+csz%i
     cend1%j = cOffset%j+csz%j
     cend1%k = cOffset%k+csz%k
     ! look over coarse grid and grab average of fine
     ii=-1
     do ic=cOffset%i+1,cend1%i
        ii=ii+2 ! 1
        jj=-1
        do jc=cOffset%j+1,cend1%j
           jj=jj+2
           kk=-1
           do kc=cOffset%k+1,cend1%k
              kk=kk+2 
#ifdef TWO_D
                 uC(ic,jc,kc,:)=.25D0*(&
                       ux(ii,  jj,kk,:)+ux(ii,  jj+1,kk,:)&
                      +ux(ii+1,jj,kk,:)+ux(ii+1,jj+1,kk,:))
                 uC2(ic,jc,kc,:)=.25D0*(&
                       ux2(ii,  jj,kk,:)+ux2(ii,  jj+1,kk,:)&
                      +ux2(ii+1,jj,kk,:)+ux2(ii+1,jj+1,kk,:))
#else 
              ! only constant restriction in 3D
              uC(ic,jc,kc,:)=.125D0*(&
                     ux(ii,  jj,kk,  :)+ux(ii,  jj+1,kk,  :)&
                    +ux(ii+1,jj,kk,  :)+ux(ii+1,jj+1,kk,  :)&
                    +ux(ii,  jj,kk+1,:)+ux(ii,  jj+1,kk+1,:)&
                    +ux(ii+1,jj,kk+1,:)+ux(ii+1,jj+1,kk+1,:))
              uC2(ic,jc,kc,:)=.125D0*(&
                    ux2(ii,  jj,kk,  :)+ux2(ii,  jj+1,kk,  :)&
                   +ux2(ii+1,jj,kk,  :)+ux2(ii+1,jj+1,kk,  :)&
                   +ux2(ii,  jj,kk+1,:)+ux2(ii,  jj+1,kk+1,:)&
                   +ux2(ii+1,jj,kk+1,:)+ux2(ii+1,jj+1,kk+1,:))
#endif
           enddo
        enddo
     enddo
  endif
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(3),ierr)
#endif
  call SetBCs(uC,cp,ct,nsg) 
  !call SetBCs(uC2,cp,ct) ! RHS does not need to be done!
  return
end subroutine RestrictFuse
!-----------------------------------------------------------------------
subroutine Prolong(uxF,uxC,ft,ct,fp,cp,cOffset,high)
  !  uxF = uxF + P * uxC
  use discretization,only:nvar,nsg
  use bc_module
  use mpistuff
  use pms,only:mg_ref_ratio
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
  !     Local vars
  integer:: ic,jc,jj,kc,ii,kk,cist,cjst,ckst,lpe
  double precision :: vv(nvar),a11,a22,a33,a44
  logical::issplit
#ifndef TWO_D
#ifdef HAVE_PETSC
  flops = 8*fp%max%hi%i*fp%max%hi%j*fp%max%hi%k
  call PetscLogEventBegin(events(5),ierr)
  call PetscLogFlops(flops,ierr)
#endif
#endif
  issplit = (fp%max%hi%i==cp%max%hi%i)
  if (mg_ref_ratio .ne. 2) stop 'Prolong: not done for refratio != 2'
  cist = cOffset%i+1
  cjst = cOffset%j+1
  ckst = cOffset%k+1
  lpe = -1
  if (issplit) then ! split grid - chop it, not SR
     if (fp%max%hi%j/=cp%max%hi%j.or.fp%max%hi%k/=cp%max%hi%k) stop 'RestrictFuse: splits not right'
     call MPI_comm_rank(ct%loc_comm,lpe,ierr)
#ifdef TWO_D
     if ( lpe == 0 ) then       ! (1,1)
     else if ( lpe == 1 ) then ! (1,2)
        cjst = cjst + cp%max%hi%j/2
     else 
        cist = cist + cp%max%hi%i/2
        if ( lpe == 2 ) then ! (2,1)
        else if ( lpe == 3 ) then ! (2,2)
           cjst = cjst + cp%max%hi%j/2
        end if
     end if
#else
     if ( lpe == 0 ) then       ! (1,1,1)
     else if ( lpe == 1 ) then ! (1,1,2)
        ckst = ckst + cp%max%hi%k/2
     else if ( lpe == 2 ) then ! (1,2,1)
        cjst = cjst + cp%max%hi%j/2
     else if ( lpe == 3 ) then ! (1,2,2)
        cjst = cjst + cp%max%hi%j/2
        ckst = ckst + cp%max%hi%k/2
     else
        cist = cist+cp%max%hi%i/2
        if ( lpe == 4 ) then  ! (2,1,1)
        else if ( lpe == 5 ) then ! (2,1,2)
           ckst = ckst + cp%max%hi%k/2
        else if ( lpe == 6 ) then ! (2,2,1)
           cjst = cjst + cp%max%hi%j/2
        else if ( lpe == 7 ) then ! (2,2,2)
           cjst = cjst + cp%max%hi%j/2
           ckst = ckst + cp%max%hi%k/2
        else
           stop 'should not be here 3): Prolong'
        endif
     end if
#endif
  endif ! else not split
#ifdef TWO_D
  a11 = 9.d0/16.d0
  a22 = 3.d0/16.d0
  a33 = 1.d0/16.d0
#else
  a11 = 27.d0/64.d0
  a22 = 9.d0/64.d0
  a33 = 3.d0/64.d0
  a44 = 1.d0/64.d0
#endif
!!$  if (.not.issplit.and.high.and.mype==4) then
!!$     write(6,'(I2,A,I2,I2,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,I2)')&
!!$          mype,'] c st:',cist,cjst,ckst,&
!!$          ', C max i:',cp%max%hi%i,', j:',cp%max%hi%j,', k:',cp%max%hi%k,&
!!$          ', F max i:',fp%max%hi%i,', j:',fp%max%hi%j,', k:',fp%max%hi%k
!!$  end if
  kc = ckst-1
  do kk=1,fp%max%hi%k,2
     kc = kc+1
     jc = cjst-1
     do jj=1,fp%max%hi%j,2
        jc = jc+1
        ic = cist-1
        do ii=1,fp%max%hi%i,2
           ic=ic+1
!!$if (.not.issplit.and.high.and.mype==4) write(6,'(I3,A,I4,I4,I4,A,I4,I4,I4,A,E15.5)')mype,&
!!$     '] fine:',ii,jj,kk,', coarse:',ic,jc,kc,', UC=',uxC(ic,jc,kc,:)
!!$     
           vv = a11*uxC(ic,jc,kc,:)
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic,jc-1,kc,:)
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           
           vv = a22*uxC(ic,jc+1,kc,:)
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic+1,jc,kc,:)
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic-1,jc,kc,:)      
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic-1,jc-1,kc,:)
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           
           vv = a33*uxC(ic-1,jc+1,kc,:)
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic+1,jc-1,kc,:)
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           
           vv = a33*uxC(ic+1,jc+1,kc,:)
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
#ifndef TWO_D
           vv = a11*uxC(ic,jc,kc,:)
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           
           vv = a22*uxC(ic,jc-1,kc,:)
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv
           
           vv = a22*uxC(ic,jc+1,kc,:)
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
           
           vv = a22*uxC(ic+1,jc,kc,:)
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
           
           vv = a22*uxC(ic-1,jc,kc,:)          
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic-1,jc-1,kc,:)
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           
           vv = a33*uxC(ic-1,jc+1,kc,:)            
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic+1,jc-1,kc,:)            
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv
           
           vv = a33*uxC(ic+1,jc+1,kc,:)            
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
                
           vv = a22*uxC(ic,jc,kc-1,:)
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic,jc,kc+1,:)
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic-1,jc,kc-1,:)
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic-1,jc,kc+1,:)
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic,jc-1,kc-1,:)
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           
           vv = a33*uxC(ic,jc-1,kc+1,:)
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv
           
           vv = a33*uxC(ic+1,jc,kc-1,:)
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic+1,jc,kc+1,:)
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic,jc+1,kc-1,:)
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic,jc+1,kc+1,:)
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
                
           vv = a44*uxC(ic-1,jc-1,kc-1,:)
           uxF(ii,jj,kk,:)=uxF(ii,jj,kk,:) + vv
           
           vv = a44*uxC(ic+1,jc-1,kc-1,:)
           uxF(ii+1,jj,kk,:)=uxF(ii+1,jj,kk,:) + vv
           
           vv = a44*uxC(ic-1,jc+1,kc-1,:)
           uxF(ii,jj+1,kk,:)=uxF(ii,jj+1,kk,:) + vv
           
           vv = a44*uxC(ic-1,jc-1,kc+1,:)
           uxF(ii,jj,kk+1,:)=uxF(ii,jj,kk+1,:) + vv
           
           vv = a44*uxC(ic+1,jc+1,kc-1,:)
           uxF(ii+1,jj+1,kk,:)=uxF(ii+1,jj+1,kk,:) + vv
           
           vv = a44*uxC(ic-1,jc+1,kc+1,:)
           uxF(ii,jj+1,kk+1,:)=uxF(ii,jj+1,kk+1,:) + vv
           
           vv = a44*uxC(ic+1,jc-1,kc+1,:)
           uxF(ii+1,jj,kk+1,:)=uxF(ii+1,jj,kk+1,:) + vv

           vv = a44*uxC(ic+1,jc+1,kc+1,:)
           uxF(ii+1,jj+1,kk+1,:)=uxF(ii+1,jj+1,kk+1,:) + vv
#endif
        enddo
     enddo
  enddo
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(5),ierr)
#endif
  call SetBCs(uxF,fp,ft,nsg)
  return
end subroutine Prolong
!-----------------------------------------------------------------------
subroutine GS_RB_const_Lap(phi,rhs,p,t,nits)
  use discretization, only:nvar,nsg,alpha
  use bc_module
  use mpistuff
  use pms, only:verbose
  implicit none
  type(patcht)::p
  double precision,dimension(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nvar)::phi,rhs
  type(topot),intent(in)::t
  integer,intent(in)::nits
  !     Local vars
  integer:: m,ii,jj,j,kk,rbi,offi,ig,jg,kg
  double precision:: ti,tj,tk,dxi2,dyi2,numer,denoi
#ifndef TWO_D
  double precision:: dzi2
#endif

  ! these are always even for non-SR and in SR does not matter except for exact
  ig = 0 ! getIglobalx(val,t) 
  jg = 0 ! getIglobaly(val,t)
  kg = 0 ! getIglobalz(val,t)
#ifndef TWO_D
  flops = 13*p%max%hi%i*p%max%hi%j*p%max%hi%k/2 
#else
  flops = 9*p%max%hi%i*p%max%hi%j*p%max%hi%k/2 
#endif
  dxi2=1.d0/p%dx%i**2
  dyi2=1.d0/p%dx%j**2
  denoi=alpha*2.d0*(dxi2 + dyi2)  ! diag
#ifndef TWO_D
  dzi2=1.d0/p%dx%k**2
  denoi=denoi + alpha*2.d0*dzi2
#endif
  denoi=-1.d0/denoi
  do m=1,nits
     ! red/black
#ifdef HAVE_PETSC
     call PetscLogEventBegin(events(2),ierr)
     call PetscLogFlops(flops,ierr)
#endif 
     do rbi = 0,1
        do kk=1,p%max%hi%k
           do jj=1,p%max%hi%j
              !offi = mod(ig+kk+jg+jj+kg-2+rbi,2)+1
              offi = mod(kk+jj+rbi,2) + 1
              do ii=offi,p%max%hi%i,2
                 ti = alpha*dxi2*(phi(ii+1,jj,kk,1)+phi(ii-1,jj,kk,1))
                 tj = alpha*dyi2*(phi(ii,jj+1,kk,1)+phi(ii,jj-1,kk,1))
                 !     set
                 numer = rhs(ii,jj,kk,1) - ti - tj
#ifndef TWO_D
                 ! Z direction
                 !     set
                 tk = alpha*dzi2*(phi(ii,jj,kk+1,1) + phi(ii,jj,kk-1,1))
                 numer = numer - tk
#endif
                 phi(ii,jj,kk,1) = numer*denoi
              enddo       ! ii
           enddo          ! jj
        enddo             ! kk
     enddo                ! r/b i
#ifdef HAVE_PETSC
     call PetscLogEventEnd(events(2),ierr)
#endif
     call SetBCs(phi,p,t,nsg)
  enddo                   ! iters
end subroutine GS_RB_const_Lap
!-----------------------------------------------------------------------
subroutine GS_Lex_const_Lap(phi,rhs,p,t,nits)
  use discretization, only:nvar,nsg,alpha
  use bc_module
  use mpistuff
  use pms, only:verbose
  implicit none
  type(patcht)::p
  double precision,dimension(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nvar)::phi,rhs
  type(topot),intent(in)::t
  integer,intent(in):: nits
  !     Local vars
  integer:: m,ii,jj,j,kk
  double precision:: ti,tj,tk,dxi2,dyi2,numer,denoi
#ifndef TWO_D
  double precision:: dzi2
#endif
#ifndef TWO_D
  flops = 13*p%max%hi%i*p%max%hi%j*p%max%hi%k/2 
#else
  flops = 9*p%max%hi%i*p%max%hi%j*p%max%hi%k/2 
#endif
  dxi2=1.d0/p%dx%i**2
  dyi2=1.d0/p%dx%j**2
  denoi=alpha*2.d0*(dxi2 + dyi2)  ! diag
#ifndef TWO_D
  dzi2=1.d0/p%dx%k**2
  denoi=denoi + alpha*2.d0*dzi2
#endif
  denoi=-1.d0/denoi
  do m=1,nits
#ifdef HAVE_PETSC
     call PetscLogEventBegin(events(2),ierr)
     call PetscLogFlops(flops,ierr)
#endif
     do kk=1,p%max%hi%k
        do jj=1,p%max%hi%j
           do ii=1,p%max%hi%i
              ti = alpha*dxi2*(phi(ii+1,jj,kk,1)+phi(ii-1,jj,kk,1))
              tj = alpha*dyi2*(phi(ii,jj+1,kk,1)+phi(ii,jj-1,kk,1))
              !     set
              numer = rhs(ii,jj,kk,1) - ti - tj                 
#ifndef TWO_D
              ! Z direction
              !     set
              tk = alpha*dzi2*(phi(ii,jj,kk+1,1)+phi(ii,jj,kk-1,1))
              numer = numer - tk
#endif
              phi(ii,jj,kk,1) = numer*denoi
           enddo       ! ii
        enddo          ! jj
     enddo             ! kk
#ifdef HAVE_PETSC
     call PetscLogEventEnd(events(2),ierr)
#endif
     call SetBCs(phi,p,t,nsg)
  end do ! its
  return
end subroutine GS_Lex_const_Lap
!-----------------------------------------------------------------------
subroutine Jacobi_const_Lap(phi,rhs,p,t,nits) 
  use discretization, only:nvar,alpha
  use base_data_module
  use mpistuff
  implicit none
  type(patcht)::p
  double precision,intent(out),dimension(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nvar)::phi
  double precision,intent(in),dimension(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nvar)::rhs
  type(topot),intent(in)::t
  integer,intent(in):: nits
  !     Local vars
  double precision,dimension(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nvar)::Res,Dk,Aux
  double precision:: diagi,dxl,dyl,ti,tj,tk,dxi2,dyi2,dzi2,dzl
  double precision:: over,under,rhok,rhokp1,beta,alpha2,delta,theta,s1,ct1,ct2,omega,norm
  integer::k
#ifndef TWO_D
  flops = 13*p%max%hi%i*p%max%hi%j*p%max%hi%k/2 
#else
  flops = 9*p%max%hi%i*p%max%hi%j*p%max%hi%k/2 
#endif
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(4),ierr)
  call PetscLogFlops(flops,ierr)
#endif
  dxl=p%dx%i
  dyl=p%dx%j
  dxi2=1.d0/dxl**2
  dyi2=1.d0/dyl**2
  diagi=-alpha*2.d0*(dxi2+dyi2)
#ifndef TWO_D
  dzl=p%dx%k
  dzi2=1.d0/dzl**2
  diagi=-alpha*2.d0*(dxi2+dyi2+dzi2)
#endif
  diagi=1.d0/diagi
  rhok = 2.d0 ! max eigen of D^-1A
  if (.true.) then
     omega = diagi*4.d0/(13.d0*rhok)
     ti=norm(phi,p%all,p%max,p%dx,t%comm,2)
if(mype==0)print *, 'Jacobi_const_Lap: omega=', omega, ', |phi|=', ti
     do k=0,nits-1
        call Apply_const_Lap(Aux,phi,p,t)
        Aux = rhs - Aux
        phi = phi + omega*Aux
ti=norm(Aux,p%all,p%max,p%dx,t%comm,2)
if(mype==0)print *, k,'|r|=',ti
     enddo
  else
     over = 1.0d0
     under = 1.d0/2.d0     
     beta =  rhok * over 
     alpha2 = rhok * under
     delta = (beta - alpha2)/2.
     theta = (beta + alpha2)/2.
     s1 = theta / delta
     rhok = 1./s1

     call Apply_const_Lap(Dk,phi,p,t)
     Dk = rhs - Dk
     Dk = Dk * diagi
     ct1 = 1.d0/theta
     Dk = Dk*ct1
     phi = phi + Dk
     
     do k=0,nits-2
        rhokp1 = 1.d0/(2.d0*s1 - rhok)
        ct1 = rhokp1*rhok
        ct2 = 2.d0*rhokp1/delta
        rhok = rhokp1

        call Apply_const_Lap(Aux,phi,p,t)
        Aux = rhs - Aux
        Res = Aux * diagi 
        
        ! Dk[] = rhokp1 * rhok * Dk[] + 2. * rhokp1 * res[] / ( delta * diag[] )
        Dk =  ct1*Dk + ct2*Res
        phi = phi + Dk
     enddo                   ! iters
  end if
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(4),ierr)
#endif
end subroutine Jacobi_const_Lap
!-----------------------------------------------------------------------
subroutine Apply_const_Lap(uxo,ux,p,t)
  use discretization, only:nvar,nsg,alpha
  use bc_module
  use mpistuff
  implicit none
  type(patcht)::p
  type(topot),intent(in)::t
  double precision,intent(out):: uxo(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nvar)
  double precision,intent(in)::ux(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nvar)
  !     Local vars
  integer::its,ii,jj,kk
  double precision::dxl,dyl,ti,tj,tk,dxi2,dyi2,dzi2,dzl
#ifndef TWO_D
  flops = 13*p%max%hi%i*p%max%hi%j*p%max%hi%k 
#else
  flops = 9*p%max%hi%i*p%max%hi%j*p%max%hi%k 
#endif
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(4),ierr)
  call PetscLogFlops(flops,ierr)
#endif
  dxl=p%dx%i
  dyl=p%dx%j
  dxi2=1.d0/dxl**2
  dyi2=1.d0/dyl**2
#ifndef TWO_D
  dzl=p%dx%k
  dzi2=1.d0/dzl**2
#endif
  do kk=1,p%max%hi%k
     do jj=1,p%max%hi%j
        do ii=1,p%max%hi%i
           ti =   dxi2*(ux(ii+1,jj,kk,1)+ux(ii-1,jj,kk,1))&
                + dyi2*(ux(ii,jj+1,kk,1)+ux(ii,jj-1,kk,1))&
                - 2.d0*(dxi2+dyi2)*ux(ii,jj,kk,1)
#ifndef TWO_D
           ti = ti + dzi2*(ux(ii,jj,kk+1,1)+ux(ii,jj,kk-1,1))&
                - 2.d0*dzi2*ux(ii,jj,kk,1)
#endif
           uxo(ii,jj,kk,1) = alpha*ti               
        enddo             
     enddo                
  enddo                   
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(4),ierr)
#endif
  call SetBCs(uxo,p,t,nsg)
  return
end subroutine Apply_const_Lap
!-----------------------------------------------------------------------
subroutine formExactU(exact,all,val,ip,dx)
  use bc_module
  use mpistuff
  use pms
  use domain
  implicit none
  type(box),intent(in)::val
  type(ipoint),intent(in)::ip
  type(box),intent(in)::all
  type(dpoint),intent(in):: dx
  double precision,intent(out)::exact(all%lo%i:all%hi%i,all%lo%j:all%hi%j,all%lo%k:all%hi%k,1)
  !
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(7),ierr)
#endif
  ig = getIglobalx(val,ip)-1
  jg = getIglobaly(val,ip)-1
  kg = getIglobalz(val,ip)-1
  ! have to do exact by hand - valid region plus ghost
  do kk=val%lo%k,val%hi%k
     do jj=val%lo%j,val%hi%j
        do ii=val%lo%i,val%hi%i
           coord(1) = xl+(ig+ii-val%lo%i)*dx%i+0.5*dx%i
           coord(2) = yl+(jg+jj-val%lo%j)*dx%j+0.5*dx%j
#ifndef TWO_D
           coord(3) = zl+(kg+kk-val%lo%k)*dx%k+0.5*dx%k
#endif
           select case (problemType)
           case(0)                            
              x2 = coord*coord
              a = x2*(x2-1.)              
#ifdef TWO_D
              exact(ii,jj,kk,1) = a(1)*a(2)
#else
              exact(ii,jj,kk,1) = a(1)*a(2)*a(3)
#endif
           case(1)
              exact = 0.d0  
              print *,'**************** formExact not done for proble type 1: todo'
           end select
        end do
     end do
  end do
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(7),ierr)
#endif
  return
end subroutine formExactU
!-----------------------------------------------------------------------
subroutine FormF(rhs,p,val,ip)
  use bc_module
  use pms
  use domain
  use mpistuff 
  use discretization, only:alpha
  implicit none
  type(patcht)::p             ! max, dx & allocated
  type(box),intent(in)::val   ! needed to get correct coordinates with buffers in max (SR)
  type(ipoint),intent(in)::ip ! global position for coordinates
  double precision,intent(out)::rhs(p%all%lo%i:p%all%hi%i,p%all%lo%j:p%all%hi%j,p%all%lo%k:p%all%hi%k,1)

  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(7),ierr)
#endif 
  ig = getIglobalx(val,ip)-1 ! zero based
  jg = getIglobaly(val,ip)-1
  kg = getIglobalz(val,ip)-1

  select case (problemType)
  case(0)
     ! x^4 - x^2
#ifdef XPERIODIC
     stop 'FormRHS: x^4 - x^2 XPERIODIC not defined'
#endif
#ifdef YPERIODIC
     stop 'FormRHS: x^4 - x^2 YPERIODIC not defined'
#endif
#ifdef ZPERIODICtod
     stop 'FormRHS: x^4 - x^2 ZPERIODIC not defined'
#endif
     ! RHS = Lap(x^4-x^2)
     do kk=1,p%max%hi%k
        do jj=1,p%max%hi%j
           do ii=1,p%max%hi%i
              coord(1) = xl+(ig+ii-val%lo%i)*p%dx%i + 0.5*p%dx%i
              coord(2) = yl+(jg+jj-val%lo%j)*p%dx%j + 0.5*p%dx%j
#ifndef TWO_D
              coord(3) = zl+(kg+kk-val%lo%k)*p%dx%k + 0.5*p%dx%k
#endif
              x2 = coord*coord
              a = x2*(x2-1.)
              b = 12.*x2-2.
#ifdef TWO_D
              rhs(ii,jj,kk,1) = b(1)*a(2) + a(1)*b(2)
#else
              rhs(ii,jj,kk,1) = b(1)*a(2)*a(3) + a(1)*b(2)*a(3) + a(1)*a(2)*b(3)
#endif
              rhs(ii,jj,kk,1) = alpha*rhs(ii,jj,kk,1)
           end do
        end do
     end do
  case(1)
     !     for bubble
#ifndef XPERIODIC
     stop 'FormRHS: bubble not XPERIODIC not defined'
#endif
#ifndef YPERIODIC
     stop 'bFormRHS: ubble not YPERIODIC defined'
#endif
#ifndef ZPERIODIC
     stop 'FormRHS: bubble not ZPERIODIC defined'
#endif
     ! RHS = bubble
     

     stop 'FormRHS not impl -- todo'

  end select
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(7),ierr)
#endif
  return
end subroutine FormF
!-----------------------------------------------------------------
! norms
!-----------------------------------------------------------------
double precision function level_norm1(ux,all,val,dx,comm)
  use base_data_module
  use mpistuff
  implicit none 
  type(box),intent(in)::val
  type(box),intent(in)::all
  type(dpoint),intent(in):: dx
  double precision,intent(in)::ux(all%lo%i:all%hi%i,all%lo%j:all%hi%j,all%lo%k:all%hi%k,1)
  integer,intent(in)::comm

  double precision::t1,t2,vol
  
#ifndef TWO_D
  vol = dx%i*dx%j*dx%k
#else
  vol = dx%i*dx%j
#endif

  t1 = sum(abs(ux(val%lo%i:val%hi%i,val%lo%j:val%hi%j,val%lo%k:val%hi%k,1)))
  call MPI_Allreduce(t1,t2,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  level_norm1 = t2*vol
end function level_norm1
!-----------------------------------------------------------------
double precision function level_norm2(ux,all,val,dx,comm)
  use base_data_module
  use mpistuff
  implicit none 
  type(box),intent(in)::val
  type(box),intent(in)::all
  type(dpoint),intent(in):: dx
  double precision,intent(in)::ux(all%lo%i:all%hi%i,all%lo%j:all%hi%j,all%lo%k:all%hi%k,1)
  integer,intent(in)::comm
  !     
  double precision::t1,t2,vol  
  integer :: i,j,k
  
#ifndef TWO_D
  vol = dx%i*dx%j*dx%k
#else
  vol = dx%i*dx%j
#endif
  t1 = sum(ux(val%lo%i:val%hi%i,val%lo%j:val%hi%j,val%lo%k:val%hi%k,1)**2)
  call MPI_Allreduce(t1,t2,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  IF (ierr /= 0) STOP "level_norm2: *** level_norm2 all reduce error ***"
  level_norm2 = sqrt(t2*vol)
end function level_norm2
!-----------------------------------------------------------------
double precision function level_norminf(ux,all,val,dx,comm)
  use base_data_module
  use mpistuff
  implicit none 
  type(box),intent(in)::val
  type(box),intent(in)::all
  type(dpoint),intent(in):: dx
  double precision,intent(in)::ux(all%lo%i:all%hi%i,all%lo%j:all%hi%j,all%lo%k:all%hi%k,1)
  integer,intent(in)::comm
  !     
  double precision::t1,t2

  t1 = maxval(abs(ux(val%lo%i:val%hi%i,val%lo%j:val%hi%j,val%lo%k:val%hi%k,1)))
  call MPI_Allreduce(t1,t2,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  level_norminf = t2
end function level_norminf
!-----------------------------------------------------------------
double precision function norm(ux,all,val,dx,comm,type)
  use base_data_module
  use mpistuff
  implicit none 
  type(box),intent(in)::val
  type(box),intent(in)::all
  type(dpoint),intent(in):: dx
  double precision,intent(in)::ux(all%lo%i:all%hi%i,all%lo%j:all%hi%j,all%lo%k:all%hi%k,1)
  integer,intent(in) :: type
  integer,intent(in)::comm
  !
  double precision:: level_norm1,level_norm2,level_norminf
  !     
  if (type==1) then
     norm = level_norm1  (ux,all,val,dx,comm)
  else if (type==2) then
     norm = level_norm2  (ux,all,val,dx,comm)
  else if (type==3) then
     norm = level_norminf(ux,all,val,dx,comm)
  else 
     stop 'norm: unknown norm type'
  end if
end function norm
