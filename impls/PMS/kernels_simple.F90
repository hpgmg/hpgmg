!-----------------------------------------------------------------
!     Ravi Samtaney & Mark Adams
!     Copyright 2014
!-----------------------------------------------------------------
double precision function level_norm1(ux,g)
  use pe_patch_data_module
  use mpistuff
  implicit none
  !     
  type(pe_patch)::g
  double precision::ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  !
  double precision::t1,t2,vol    
  
  vol = g%dxg*g%dyg &
#ifndef TWO_D
       *g%dzg
#endif
  t1 = sum(abs(ux(g%ivallo:g%ivalhi,g%jvallo:g%jvalhi,g%kvallo:g%kvalhi,1)))
  call MPI_Allreduce(t1,t2,1,MPI_DOUBLE_PRECISION,MPI_SUM,g%comm3D,ierr)
  level_norm1 = t2*vol
end function level_norm1
!-----------------------------------------------------------------
double precision function level_norm2(ux,g)
  use pe_patch_data_module
  use mpistuff
  implicit none
  !     
  type(pe_patch)::g
  double precision::ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  !     
  double precision::t1,t2,vol  
  integer :: i,j,k

  vol = g%dxg*g%dyg&
#ifndef TWO_D
       *g%dzg
#endif

  t1 = sum(ux(g%ivallo:g%ivalhi,g%jvallo:g%jvalhi,g%kvallo:g%kvalhi,1)**2)
  call MPI_Allreduce(t1,t2,1,MPI_DOUBLE_PRECISION,MPI_SUM,g%comm3D,ierr)
  level_norm2 = sqrt(t2*vol)

end function level_norm2
!-----------------------------------------------------------------
double precision function level_norminf(ux,g)
  use pe_patch_data_module
  use mpistuff
  implicit none
  !     
  type(pe_patch)::g
  double precision::ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  !     
  double precision::t1,t2

  t1 = maxval(abs(ux(g%ivallo:g%ivalhi,g%jvallo:g%jvalhi,g%kvallo:g%kvalhi,1)))
  call MPI_Allreduce(t1,t2,1,MPI_DOUBLE_PRECISION,MPI_MAX,g%comm3D,ierr)
  level_norminf = t2
end function level_norminf
!-----------------------------------------------------------------
double precision function norm(ux,g,type)
  !
  use pe_patch_data_module
  implicit none
  double precision level_norm2,level_norminf,level_norm1
  !     
  type(pe_patch)::g
  double precision::ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  integer :: type
  !     
  if (type==1) then
     norm = level_norm1(ux,g)
  else if (type==2) then
     norm = level_norm2(ux,g)
  else if (type==3) then
     norm = level_norminf(ux,g)
  else 
     stop 'unknown norm type'
  end if
end function norm
!-----------------------------------------------------------------------
subroutine Restrict(gc,gf,uxC,ux)
  use GridModule
  use mpistuff
  use tags
  use pms, only:mg_ref_ratio
  implicit none
  type(pe_patch):: gf, gc
  double precision,intent(out):: uxC(gc%ilo:gc%ihi, gc%jlo:gc%jhi,& 
       gc%klo:gc%khi, nvar)
  double precision,intent(in):: ux(gf%ilo:gf%ihi, gf%jlo:gf%jhi,&
       gf%klo:gf%khi, nvar)
  !integer,intent(in) :: order ! not used in 3D
  !     Local vars
#ifdef TWO_D
  double precision,dimension(gc%imax/2, gc%jmax/2, 1,         nvar)::&
       sb,rb1,rb2,rb3,rb0
#else
  double precision,dimension(gc%imax/2, gc%jmax/2, gc%kmax/2, nvar)::&
       sb,rb1,rb2,rb3,rb0,rb4,rb5,rb6,rb7
#endif
  double precision :: a11,a22,a33,a44
  integer:: ic,jc,kc,cimax2,ckmax2,cjmax2,pp,lnp
  integer:: msg_id_recv(0:7),bufsz,msg_id_send(0:7)
  integer:: ii,jj,kk,ist,jst,kst
  
  if (mg_ref_ratio .ne. 2) stop 'Restrict not done for refratio != 2'

#ifdef HAVE_PETSC
  flops = 2*gf%imax*gf%jmax*gf%kmax
  call PetscLogEventBegin(events(3),ierr)
  call PetscLogFlops(flops,ierr)
#endif
  
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
  
  if(gf%imax == gc%imax)then ! split grid
     ist = gc%imax/2+1
     jst = gc%jmax/2+1
#ifdef TWO_D
     bufsz = (nvar)*(gc%jmax/2)*(gc%imax/2)
     lnp = 4
     ckmax2 = 1
#else 
     bufsz = (nvar)*(gc%jmax/2)*(gc%imax/2)*(gc%kmax/2)
     lnp = 8
     ckmax2 = gc%kmax/2
#endif
     call MPI_Irecv(rb0,bufsz,MPI_DOUBLE_PRECISION,&
          0,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(0),ierr)
     call MPI_Irecv(rb1,bufsz,MPI_DOUBLE_PRECISION,&
          1,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(1),ierr)
     call MPI_Irecv(rb2,bufsz,MPI_DOUBLE_PRECISION,&
          2,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(2),ierr)
     call MPI_Irecv(rb3,bufsz,MPI_DOUBLE_PRECISION,&
          3,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(3),ierr)
#ifndef TWO_D  
     call MPI_Irecv(rb4,bufsz,MPI_DOUBLE_PRECISION,&
          4,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(4),ierr)
     call MPI_Irecv(rb5,bufsz,MPI_DOUBLE_PRECISION,&
          5,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(5),ierr)
     call MPI_Irecv(rb6,bufsz,MPI_DOUBLE_PRECISION,&
          6,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(6),ierr)
     call MPI_Irecv(rb7,bufsz,MPI_DOUBLE_PRECISION,&
          7,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(7),ierr)
#endif
     !     only one quarter/eighth of coarse vector (all of fine)
     cimax2 = gc%imax/2
     cjmax2 = gc%jmax/2
     do kc=1,ckmax2,1
        kk=(kc-1)*2+1
        do jc=1,cjmax2,1
           jj=(jc-1)*2+1
           do ic=1,cimax2,1
              ii=(ic-1)*2+1
#ifdef TWO_D
              if(.true.) then
                 sb(ic,jc,kc,:)=.25d0*(&
                      ux(ii,jj,kk,:)+ux(ii+1,jj,kk,:)&
                      +ux(ii+1,jj+1,kk,:)+ux(ii,jj+1,kk,:))
              else
                 sb(ic,jc,kc,:)=a11*(&
                      ux(ii,jj,kk,:)+ux(ii,jj+1,kk,:)&
                      +ux(ii+1,jj,kk,:)+ux(ii+1,jj+1,kk,:)) +&
                      a22*(&
                      ux(ii-1,jj,kk,:) +    ux(ii+2,jj,kk,:)&
                      +ux(ii-1,jj+1,kk,:) + ux(ii+2,jj+1,kk,:)&
                      +ux(ii,jj-1,kk,:) +   ux(ii,jj+2,kk,:)&
                      +ux(ii+1,jj-1,kk,:) + ux(ii+1,jj+2,kk,:)) +&
                      a33*(&
                      ux(ii-1,jj-1,kk,:)  + ux(ii+2,jj-1,kk,:)&
                      +ux(ii-1,jj+2,kk,:) + ux(ii+2,jj+2,kk,:))
              endif
#else 
              if(.true.) then
                 sb(ic,jc,kc,:)=.125d0*(&
                       ux(ii,  jj,  kk,  :)+ux(ii+1,jj,  kk,  :)&
                      +ux(ii+1,jj+1,kk,  :)+ux(ii,  jj+1,kk,  :)&
                      +ux(ii,  jj,  kk+1,:)+ux(ii+1,jj,  kk+1,:)&
                      +ux(ii+1,jj+1,kk+1,:)+ux(ii,  jj+1,kk+1,:))
              else
                 stop '3D linear not done Restrict'
              endif
#endif
           enddo
        enddo
     enddo
     ! send to my cohort & me
     do pp=0,lnp-1            
        call MPI_Isend(sb,bufsz,MPI_DOUBLE_PRECISION,pp,&
             MSG_RESTRICT_TAG,gc%loc_comm,msg_id_send(pp),ierr)
     enddo
     MSG_RESTRICT_TAG = MSG_RESTRICT_TAG + 1
     !     0    recv (includes self recv)  (1,1[,1])
     call MPI_Wait(msg_id_recv(0), status, ierr)
     do kc=1,ckmax2,1
        do jc=1,cjmax2,1
           do ic=1,cimax2,1
              uxC(ic,jc,kc,:) = rb0(ic,jc,kc,:)
           enddo
        enddo
     enddo
     !     1
     call MPI_Wait(msg_id_recv(1), status, ierr)
#ifdef TWO_D    
     kc = 1; jj = 1         ! (1,2)
     do jc=jst,gc%jmax,1
        do ic=1,cimax2,1
           uxC(ic,jc,kc,:) = rb1(ic,jj,kc,:)
        enddo
        jj = jj + 1
     enddo
#else 
     kst = gc%kmax/2+1      
     kk = 1
     do kc=kst,gc%kmax,1    ! (1,1,2)
        do jc=1,cjmax2,1
           do ic=1,cimax2,1
              uxC(ic,jc,kc,:) = rb1(ic,jc,kk,:)
           enddo
        enddo
        kk = kk + 1
     enddo
#endif
     !     2  
     call MPI_Wait(msg_id_recv(2), status, ierr)
#ifdef TWO_D
     do jc=1,cjmax2,1       ! (2,1)
        ii = 1
        do ic=ist,gc%imax,1
           uxC(ic,jc,kc,:) = rb2(ii,jc,kc,:)
           ii = ii + 1
        enddo
     enddo
#else 
     do kc=1,ckmax2,1       ! (1,2,1)
        jj = 1 
        do jc=jst,gc%jmax,1
           do ic=1,cimax2,1
              uxC(ic,jc,kc,:) = rb2(ic,jj,kc,:)
           enddo
           jj = jj + 1
        enddo
     enddo
#endif
     !     3
     call MPI_Wait(msg_id_recv(3), status, ierr)
#ifdef TWO_D
     jj = 1                 ! (2,2)
     do jc=jst,gc%jmax,1
        ii = 1
        do ic=ist,gc%imax,1
           uxC(ic,jc,kc,:) = rb3(ii,jj,kc,:)
           ii = ii + 1
        enddo
        jj = jj + 1
     enddo
#else
     kk = 1
     do kc=kst,gc%kmax,1    ! (1,2,2)
        jj = 1
        do jc=jst,gc%jmax,1
           do ic=1,cimax2,1
              uxC(ic,jc,kc,:) = rb3(ic,jj,kk,:)
           enddo
           jj = jj + 1
        enddo
        kk = kk + 1
     enddo
     !     4 (2,1,1)
     call MPI_Wait(msg_id_recv(4), status, ierr)
     do kc=1,ckmax2,1
        do jc=1,cjmax2,1
           ii = 1
           do ic=ist,gc%imax,1
              uxC(ic,jc,kc,:) = rb4(ii,jc,kc,:)
              ii = ii + 1
           enddo
        enddo
        
     enddo
     !     5 (2,1,2)
     call MPI_Wait(msg_id_recv(5), status, ierr)
     kk = 1
     do kc=kst,gc%kmax,1
        do jc=1,cjmax2,1
           ii = 1
           do ic=ist,gc%imax,1
              uxC(ic,jc,kc,:) = rb5(ii,jc,kk,:)
              ii = ii + 1
           enddo
        enddo
        kk = kk + 1
     enddo
     !     6 (2,2,1)
     call MPI_Wait(msg_id_recv(6), status, ierr)
     do kc=1,ckmax2,1
        jj = 1
        do jc=jst,gc%jmax,1
           ii = 1
           do ic=ist,gc%imax,1
              uxC(ic,jc,kc,:) = rb6(ii,jj,kc,:)
              ii = ii + 1
           enddo
           jj = jj + 1
        enddo
     enddo
     !     7 (2,2,2)
     call MPI_Wait(msg_id_recv(7), status, ierr)
     kk = 1
     do kc=kst,gc%kmax,1
        jj = 1
        do jc=jst,gc%jmax,1
           ii = 1
           do ic=ist,gc%imax,1
              uxC(ic,jc,kc,:) = rb7(ii,jj,kk,:)
              ii = ii + 1
           enddo
           jj = jj + 1
        enddo
        kk = kk + 1
     enddo
#endif     
     do pp=0,lnp-1
        call MPI_Wait(msg_id_send(pp), status, ierr)
     enddo
  else ! normal grid (not split)
     do kc=1,gc%kmax,1
        kk=(kc-1)*2+1 
        do jc=1,gc%jmax,1
           jj=(jc-1)*2+1
           do ic=1,gc%imax,1
              ii=(ic-1)*2+1
#ifdef TWO_D
              if(.true.)then
                 uxC(ic,jc,kc,:)=.25D0*(&
                      ux(ii,jj,kc,:)+ux(ii,jj+1,kc,:)&
                      +ux(ii+1,jj,kc,:)+ux(ii+1,jj+1,kc,:))
              else
                 uxC(ic,jc,kc,:)=a11*(&  
                      ux(ii,jj,kc,:)+ux(ii,jj+1,kc,:)&
                      +ux(ii+1,jj,kc,:)+ux(ii+1,jj+1,kc,:)) +&
                      a22*(&
                      ux(ii-1,jj,kc,:)+ux(ii+2,jj,kc,:)&
                      +ux(ii-1,jj+1,kc,:)+ux(ii+2,jj+1,kc,:)) +&
                      a22*(&
                      ux(ii,jj-1,kc,:)+ux(ii,jj+2,kc,:)&
                      +ux(ii+1,jj-1,kc,:)+ux(ii+1,jj+2,kc,:)) +&
                      a33*(&
                      ux(ii-1,jj-1,kc,:)+ux(ii+2,jj-1,kc,:)&
                      +ux(ii-1,jj+2,kc,:)+ux(ii+2,jj+2,kc,:))
              endif
#else 
              ! only constant restriction in 3D
              uxC(ic,jc,kc,:)=.125D0*(&
                    ux(ii,  jj,kk,  :)+ux(ii,  jj+1,kk,  :)&
                   +ux(ii+1,jj,kk,  :)+ux(ii+1,jj+1,kk,  :)&
                   +ux(ii,  jj,kk+1,:)+ux(ii,  jj+1,kk+1,:)&
                   +ux(ii+1,jj,kk+1,:)+ux(ii+1,jj+1,kk+1,:))
#endif
           enddo
        enddo
     enddo
  endif
    
  call SetBCs(uxC,gc)

!!$  call flush(6)
!!$  call sleep(mype)
!!$  do kk=1,gc%kmax
!!$     write(6,*) '[',mype,'] Restrict done:k=',kk
!!$     do jj=1,gc%jmax
!!$        write(6,'(E11.3,E11.3,E11.3,E11.3)') uxC(1:gc%imax,jj,kk,1)
!!$     end do
!!$  end do
!!$  call flush(6)
!!$  call sleep(8-mype)
!!$  stop
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(3),ierr)
#endif

  return
end subroutine Restrict
!-----------------------------------------------------------------------
! Restrict2: fuse two restricts.  Only do BCs on 2nd!
!-----------------------------------------------------------------------
subroutine RestrictFuse(gc,gf,uC,uC2,ux,ux2)
  use GridModule
  use mpistuff
  use tags
  use pms, only:mg_ref_ratio
  implicit none
  type(pe_patch):: gf,gc
  double precision,intent(in):: ux(gf%ilo:gf%ihi, gf%jlo:gf%jhi,&
       gf%klo:gf%khi, nvar)
  double precision,intent(in):: ux2(gf%ilo:gf%ihi, gf%jlo:gf%jhi,&
       gf%klo:gf%khi, nvar)
  double precision,intent(out):: uC(gc%ilo:gc%ihi, gc%jlo:gc%jhi,& 
       gc%klo:gc%khi, nvar)
  double precision,intent(out):: uC2(gc%ilo:gc%ihi, gc%jlo:gc%jhi,& 
       gc%klo:gc%khi, nvar)
  !integer,intent(in) :: order ! not used in 3D
  !     Local vars
#ifdef TWO_D
  double precision,dimension(gc%imax/2,gc%jmax/2,1,2*nvar)::&
       sb,rb1,rb2,rb3,rb0
#else
  double precision,dimension(gc%imax/2,gc%jmax/2,gc%kmax/2,2*nvar)::&
       sb,rb0,rb1,rb2,rb3,rb4,rb5,rb6,rb7
#endif
  double precision :: a11,a22,a33,a44,tmp
  integer:: ic,jc,kc,cimax2,ckmax2,cjmax2,pp,lnp
  integer:: msg_id_recv(0:7),bufsz,msg_id_send(0:7)
  integer:: ii,jj,kk,ist,jst,kst
  
  if (mg_ref_ratio .ne. 2) stop 'Restrict2 not done for refratio != 2'

#ifdef HAVE_PETSC
  flops = 2*2*gf%imax*gf%jmax*gf%kmax
  call PetscLogEventBegin(events(3),ierr)
  call PetscLogFlops(flops,ierr)
#endif

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
  if(gf%imax == gc%imax)then ! split grid
     ist = gc%imax/2+1
     jst = gc%jmax/2+1
#ifdef TWO_D
     bufsz = 2*nvar*(gc%jmax/2)*(gc%imax/2)
     lnp = 4
     ckmax2 = 1
#else 
     bufsz = 2*nvar*(gc%jmax/2)*(gc%imax/2)*(gc%kmax/2)
     lnp = 8
     ckmax2 = gc%kmax/2
#endif
     call MPI_Irecv(rb0,bufsz,MPI_DOUBLE_PRECISION,&
          0,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(0),ierr)
     call MPI_Irecv(rb1,bufsz,MPI_DOUBLE_PRECISION,&
          1,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(1),ierr)
     call MPI_Irecv(rb2,bufsz,MPI_DOUBLE_PRECISION,&
          2,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(2),ierr)
     call MPI_Irecv(rb3,bufsz,MPI_DOUBLE_PRECISION,&
          3,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(3),ierr)
#ifndef TWO_D  
     call MPI_Irecv(rb4,bufsz,MPI_DOUBLE_PRECISION,&
          4,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(4),ierr)
     call MPI_Irecv(rb5,bufsz,MPI_DOUBLE_PRECISION,&
          5,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(5),ierr)
     call MPI_Irecv(rb6,bufsz,MPI_DOUBLE_PRECISION,&
          6,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(6),ierr)
     call MPI_Irecv(rb7,bufsz,MPI_DOUBLE_PRECISION,&
          7,MSG_RESTRICT_TAG,gc%loc_comm,msg_id_recv(7),ierr)
#endif
     !     only one quarter/eighth of coarse vector (all of fine)
     cimax2 = gc%imax/2
     cjmax2 = gc%jmax/2
     do kc=1,ckmax2,1
        kk=(kc-1)*2+1
        do jc=1,cjmax2,1
           jj=(jc-1)*2+1
           do ic=1,cimax2,1
              ii=(ic-1)*2+1
#ifdef TWO_D
              if(.true.) then
                 sb(ic,jc,kc,1:nvar)=.25d0*(&
                      ux(ii,jj,kk,:)+ux(ii+1,jj,kk,:)&
                      +ux(ii+1,jj+1,kk,:)+ux(ii,jj+1,kk,:))
                 sb(ic,jc,kc,nvar+1:2*nvar)=.25d0*(&
                      ux2(ii,jj,kk,:)+ux2(ii+1,jj,kk,:)&
                      +ux2(ii+1,jj+1,kk,:)+ux2(ii,jj+1,kk,:))
              else
                 sb(ic,jc,kc,1:nvar)=a11*(&
                      ux(ii,jj,kk,:)+ux(ii,jj+1,kk,:)&
                      +ux(ii+1,jj,kk,:)+ux(ii+1,jj+1,kk,:)) +&
                      a22*(&
                      ux(ii-1,jj,kk,:) +    ux(ii+2,jj,kk,:)&
                      +ux(ii-1,jj+1,kk,:) + ux(ii+2,jj+1,kk,:)&
                      +ux(ii,jj-1,kk,:) +   ux(ii,jj+2,kk,:)&
                      +ux(ii+1,jj-1,kk,:) + ux(ii+1,jj+2,kk,:)) +&
                      a33*(&
                      ux(ii-1,jj-1,kk,:)  + ux(ii+2,jj-1,kk,:)&
                      +ux(ii-1,jj+2,kk,:) + ux(ii+2,jj+2,kk,:))
                 sb(ic,jc,kc,nvar+1:2*nvar)=a11*(&
                      ux2(ii,jj,kk,:)+ux2(ii,jj+1,kk,:)&
                      +ux2(ii+1,jj,kk,:)+ux2(ii+1,jj+1,kk,:)) +&
                      a22*(&
                      ux2(ii-1,jj,kk,:) +    ux2(ii+2,jj,kk,:)&
                      +ux2(ii-1,jj+1,kk,:) + ux2(ii+2,jj+1,kk,:)&
                      +ux2(ii,jj-1,kk,:) +   ux2(ii,jj+2,kk,:)&
                      +ux2(ii+1,jj-1,kk,:) + ux2(ii+1,jj+2,kk,:)) +&
                      a33*(&
                      ux2(ii-1,jj-1,kk,:)  + ux2(ii+2,jj-1,kk,:)&
                      +ux2(ii-1,jj+2,kk,:) + ux2(ii+2,jj+2,kk,:))
              endif
#else 
              if(.true.) then
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
              else
                 stop '3D linear not done Restrict'
              endif
#endif
           enddo
        enddo
     enddo
     ! send to my cohort & me
     do pp=0,lnp-1            
        call MPI_Isend(sb,bufsz,MPI_DOUBLE_PRECISION,pp,&
             MSG_RESTRICT_TAG,gc%loc_comm,msg_id_send(pp),ierr)
     enddo
     MSG_RESTRICT_TAG = MSG_RESTRICT_TAG + 1
     !     0    recv (includes self recv)  (1,1[,1])
     call MPI_Wait(msg_id_recv(0), status, ierr)
     do kc=1,ckmax2,1
        do jc=1,cjmax2,1
           do ic=1,cimax2,1
              uC(ic,jc,kc,:) = rb0(ic,jc,kc,1:nvar)
              uC2(ic,jc,kc,:) = rb0(ic,jc,kc,nvar+1:2*nvar)
           enddo
        enddo
     enddo
     !     1
     call MPI_Wait(msg_id_recv(1), status, ierr)
#ifdef TWO_D    
     kc = 1; jj = 1         ! (1,2)
     do jc=jst,gc%jmax,1
        do ic=1,cimax2,1
           uC(ic,jc,kc,:) = rb1(ic,jj,kc,1:nvar)
           uC2(ic,jc,kc,:) = rb1(ic,jj,kc,nvar+1:2*nvar)
        enddo
        jj = jj + 1
     enddo
#else 
     kst = gc%kmax/2+1      
     kk = 1
     do kc=kst,gc%kmax,1    ! (1,1,2)
        do jc=1,cjmax2,1
           do ic=1,cimax2,1
              uC(ic,jc,kc,:) = rb1(ic,jc,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb1(ic,jc,kk,nvar+1:2*nvar)
           enddo
        enddo
        kk = kk + 1
     enddo
#endif
     !     2  
     call MPI_Wait(msg_id_recv(2), status, ierr)
#ifdef TWO_D
     do jc=1,cjmax2,1       ! (2,1)
        ii = 1
        do ic=ist,gc%imax,1
           uC(ic,jc,kc,:) = rb2(ii,jc,kc,1:nvar)
           uC2(ic,jc,kc,:) = rb2(ii,jc,kc,nvar+1:2*nvar)
           ii = ii + 1
        enddo
     enddo
#else 
     do kc=1,ckmax2,1       ! (1,2,1)
        jj = 1 
        do jc=jst,gc%jmax,1
           do ic=1,cimax2,1
              uC(ic,jc,kc,:) = rb2(ic,jj,kc,1:nvar)
              uC2(ic,jc,kc,:) = rb2(ic,jj,kc,nvar+1:2*nvar)
           enddo
           jj = jj + 1
        enddo
     enddo
#endif
     !     3
     call MPI_Wait(msg_id_recv(3), status, ierr)
#ifdef TWO_D
     jj = 1                 ! (2,2)
     do jc=jst,gc%jmax,1
        ii = 1
        do ic=ist,gc%imax,1
           uC(ic,jc,kc,:) = rb3(ii,jj,kc,1:nvar)
           uC2(ic,jc,kc,:) = rb3(ii,jj,kc,nvar+1:2*nvar)
           ii = ii + 1
        enddo
        jj = jj + 1
     enddo
#else
     kk = 1
     do kc=kst,gc%kmax,1    ! (1,2,2)
        jj = 1
        do jc=jst,gc%jmax,1
           do ic=1,cimax2,1
              uC(ic,jc,kc,:) = rb3(ic,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb3(ic,jj,kk,nvar+1:2*nvar)
           enddo
           jj = jj + 1
        enddo
        kk = kk + 1
     enddo
     !     4 (2,1,1)
     call MPI_Wait(msg_id_recv(4), status, ierr)
     do kc=1,ckmax2,1
        do jc=1,cjmax2,1
           ii = 1
           do ic=ist,gc%imax,1
              uC(ic,jc,kc,:) = rb4(ii,jc,kc,1:nvar)
              uC2(ic,jc,kc,:) = rb4(ii,jc,kc,nvar+1:2*nvar)
              ii = ii + 1
           enddo
        enddo
        
     enddo
     !     5 (2,1,2)
     call MPI_Wait(msg_id_recv(5), status, ierr)
     kk = 1
     do kc=kst,gc%kmax,1
        do jc=1,cjmax2,1
           ii = 1
           do ic=ist,gc%imax,1
              uC(ic,jc,kc,:) = rb5(ii,jc,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb5(ii,jc,kk,nvar+1:2*nvar)
              ii = ii + 1
           enddo
        enddo
        kk = kk + 1
     enddo
     !     6 (2,2,1)
     call MPI_Wait(msg_id_recv(6), status, ierr)
     do kc=1,ckmax2,1
        jj = 1
        do jc=jst,gc%jmax,1
           ii = 1
           do ic=ist,gc%imax,1
              uC(ic,jc,kc,:) = rb6(ii,jj,kc,1:nvar)
              uC2(ic,jc,kc,:) = rb6(ii,jj,kc,nvar+1:2*nvar)
              ii = ii + 1
           enddo
           jj = jj + 1
        enddo
     enddo
     !     7 (2,2,2)
     call MPI_Wait(msg_id_recv(7), status, ierr)
     kk = 1
     do kc=kst,gc%kmax,1
        jj = 1
        do jc=jst,gc%jmax,1
           ii = 1
           do ic=ist,gc%imax,1
              uC(ic,jc,kc,:) = rb7(ii,jj,kk,1:nvar)
              uC2(ic,jc,kc,:) = rb7(ii,jj,kk,nvar+1:2*nvar)
              ii = ii + 1
           enddo
           jj = jj + 1
        enddo
        kk = kk + 1
     enddo
#endif     
     do pp=0,lnp-1
        call MPI_Wait(msg_id_send(pp), status, ierr)
     enddo
  else ! normal grid (not split)
     do kc=1,gc%kmax,1
        kk=(kc-1)*2+1 
        do jc=1,gc%jmax,1
           jj=(jc-1)*2+1
           do ic=1,gc%imax,1
              ii=(ic-1)*2+1
#ifdef TWO_D
              if(.true.)then
                 uC(ic,jc,kc,:)=.25D0*(&
                      ux(ii,jj,kc,:)+ux(ii,jj+1,kc,:)&
                      +ux(ii+1,jj,kc,:)+ux(ii+1,jj+1,kc,:))
                 uC2(ic,jc,kc,:)=.25D0*(&
                      ux2(ii,jj,kc,:)+ux2(ii,jj+1,kc,:)&
                      +ux2(ii+1,jj,kc,:)+ux2(ii+1,jj+1,kc,:))
              else
                 uC(ic,jc,kc,:)=a11*(&  
                      ux(ii,jj,kc,:)+ux(ii,jj+1,kc,:)&
                      +ux(ii+1,jj,kc,:)+ux(ii+1,jj+1,kc,:)) +&
                      a22*(&
                      ux(ii-1,jj,kc,:)+ux(ii+2,jj,kc,:)&
                      +ux(ii-1,jj+1,kc,:)+ux(ii+2,jj+1,kc,:)) +&
                      a22*(&
                      ux(ii,jj-1,kc,:)+ux(ii,jj+2,kc,:)&
                      +ux(ii+1,jj-1,kc,:)+ux(ii+1,jj+2,kc,:)) +&
                      a33*(&
                      ux(ii-1,jj-1,kc,:)+ux(ii+2,jj-1,kc,:)&
                      +ux(ii-1,jj+2,kc,:)+ux(ii+2,jj+2,kc,:))
                 uC2(ic,jc,kc,:)=a11*(&  
                      ux2(ii,jj,kc,:)+ux2(ii,jj+1,kc,:)&
                      +ux2(ii+1,jj,kc,:)+ux2(ii+1,jj+1,kc,:)) +&
                      a22*(&
                      ux2(ii-1,jj,kc,:)+ux2(ii+2,jj,kc,:)&
                      +ux2(ii-1,jj+1,kc,:)+ux2(ii+2,jj+1,kc,:)) +&
                      a22*(&
                      ux2(ii,jj-1,kc,:)+ux2(ii,jj+2,kc,:)&
                      +ux2(ii+1,jj-1,kc,:)+ux2(ii+1,jj+2,kc,:)) +&
                      a33*(&
                      ux2(ii-1,jj-1,kc,:)+ux2(ii+2,jj-1,kc,:)&
                      +ux2(ii-1,jj+2,kc,:)+ux2(ii+2,jj+2,kc,:))
              endif
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
  !call SetBCs(uC,gc) ! RHS!
  call SetBCs(uC2,gc)
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(3),ierr)
#endif
  return
end subroutine RestrictFuse
!-----------------------------------------------------------------------
subroutine Prolong_2(ux,uxC,gf,gc)
  !  ux = ux + P * uxC
  use GridModule
  use mpistuff
  use pms, only:mg_ref_ratio
  implicit none
  type(pe_patch):: gf, gc
  double precision,intent(in):: &
       uxC(gc%ilo:gc%ihi, gc%jlo:gc%jhi, &
       gc%klo:gc%khi, nvar)
  double precision:: &
       ux(gf%ilo:gf%ihi, gf%jlo:gf%jhi, &
       gf%klo:gf%khi, nvar)
  !     Local vars
  integer:: ic,jc,jj,kc,ii,kk,ist,iend,jst,jend,kst,kend,lpe
  double precision :: vv(nvar),a11,a22,a33,a44
#ifndef TWO_D
#ifdef HAVE_PETSC
  flops = 106*gc%imax*gc%jmax*gc%kmax
  call PetscLogEventBegin(events(5),ierr)
  call PetscLogFlops(flops,ierr)
#endif
#endif
  
  if (mg_ref_ratio .ne. 2) stop 'Prolong_2 not done for refratio != 2'

  ! the whole coarse path
  ist = 1
  iend = gc%imax
  jst = 1 
  jend = gc%jmax
  kst = 1
  kend = gc%kmax
  if(gf%imax == gc%imax) then ! split grid - chop it
     call MPI_comm_rank(gc%loc_comm,lpe,ierr)
#ifdef TWO_D
     kend = 1
     if( lpe == 0 ) then       ! (1,1)
        iend = iend/2
        jend = jend/2
     else if ( lpe == 1 ) then ! (1,2)
        iend = iend/2
        jst = gc%jmax/2+1
     else if ( lpe == 2 ) then ! (2,1)
        ist = gc%imax/2+1
        jend = jend/2
     else if ( lpe == 3 ) then ! (2,2)
        ist = gc%imax/2+1
        jst = gc%jmax/2+1
     end if
#else
     if( lpe == 0 ) then       ! (1,1,1)
        iend = iend/2
        jend = jend/2
        kend = kend/2
     else if ( lpe == 1 ) then ! (1,1,2)
        iend = iend/2
        jend = jend/2
        kst = gc%kmax/2+1
     else if ( lpe == 2 ) then ! (1,2,1)
        iend = iend/2
        jst = gc%jmax/2+1
        kend = kend/2
     else if ( lpe == 3 ) then ! (1,2,2)
        iend = iend/2
        jst = gc%jmax/2+1
        kst = gc%kmax/2+1
     else
        ist = gc%imax/2+1
        if( lpe == 4 ) then  ! (2,1,1)
           jend = jend/2
           kend = kend/2
        else if ( lpe == 5 ) then ! (2,1,2)
           jend = jend/2
           kst = gc%kmax/2+1
        else if ( lpe == 6 ) then ! (2,2,1)
           jst = gc%jmax/2+1
           kend = kend/2
        else if ( lpe == 7 ) then ! (2,2,2)
           jst = gc%jmax/2+1
           kst = gc%kmax/2+1
        else
           stop 'should not be here 3): Prolong_2'
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
  do kk = 1,gf%kmax,2
     kc = kst + (kk-1)/2
     do jj = 1,gf%jmax,2
        jc = jst + (jj-1)/2
        do ii = 1,gf%imax,2
           ic = ist + (ii-1)/2
           
           vv = a11*uxC(ic,jc,kc,:);
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic,jc-1,kc,:);
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           
           vv = a22*uxC(ic,jc+1,kc,:);
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic+1,jc,kc,:);
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic-1,jc,kc,:);            
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic-1,jc-1,kc,:);
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           
           vv = a33*uxC(ic-1,jc+1,kc,:);            
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic+1,jc-1,kc,:);            
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           
           vv = a33*uxC(ic+1,jc+1,kc,:);            
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
#ifndef TWO_D
           vv = a11*uxC(ic,jc,kc,:);
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           
           vv = a22*uxC(ic,jc-1,kc,:);
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           
           vv = a22*uxC(ic,jc+1,kc,:);
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
           
           vv = a22*uxC(ic+1,jc,kc,:);
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
           
           vv = a22*uxC(ic-1,jc,kc,:);            
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic-1,jc-1,kc,:);
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           
           vv = a33*uxC(ic-1,jc+1,kc,:);            
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic+1,jc-1,kc,:);            
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           
           vv = a33*uxC(ic+1,jc+1,kc,:);            
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
                
           vv = a22*uxC(ic,jc,kc-1,:);
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
           
           vv = a22*uxC(ic,jc,kc+1,:);
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic-1,jc,kc-1,:);
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic-1,jc,kc+1,:);
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic,jc-1,kc-1,:);
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           
           vv = a33*uxC(ic,jc-1,kc+1,:);
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           
           vv = a33*uxC(ic+1,jc,kc-1,:);
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic+1,jc,kc+1,:);
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
           
           vv = a33*uxC(ic,jc+1,kc-1,:);
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
           
           vv = a33*uxC(ic,jc+1,kc+1,:);
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
                
           vv = a44*uxC(ic-1,jc-1,kc-1,:);
           ux(ii,jj,kk,:)=ux(ii,jj,kk,:) + vv
           
           vv = a44*uxC(ic+1,jc-1,kc-1,:);
           ux(ii+1,jj,kk,:)=ux(ii+1,jj,kk,:) + vv
           
           vv = a44*uxC(ic-1,jc+1,kc-1,:);
           ux(ii,jj+1,kk,:)=ux(ii,jj+1,kk,:) + vv
           
           vv = a44*uxC(ic-1,jc-1,kc+1,:);
           ux(ii,jj,kk+1,:)=ux(ii,jj,kk+1,:) + vv
           
           vv = a44*uxC(ic+1,jc+1,kc-1,:);
           ux(ii+1,jj+1,kk,:)=ux(ii+1,jj+1,kk,:) + vv
           
           vv = a44*uxC(ic-1,jc+1,kc+1,:);
           ux(ii,jj+1,kk+1,:)=ux(ii,jj+1,kk+1,:) + vv
           
           vv = a44*uxC(ic+1,jc-1,kc+1,:);
           ux(ii+1,jj,kk+1,:)=ux(ii+1,jj,kk+1,:) + vv
           
           vv = a44*uxC(ic+1,jc+1,kc+1,:);
           ux(ii+1,jj+1,kk+1,:)=ux(ii+1,jj+1,kk+1,:) + vv
#endif
        enddo
     enddo
  enddo
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(5),ierr)
#endif
  call SetBCs(ux,gf)

  return
end subroutine Prolong_2

!-----------------------------------------------------------------------
subroutine GSRB_const_Lap(phi,rhs,g,nits)
  use GridModule
  use mpistuff
  implicit none
  type(pe_patch):: g
  double precision:: phi(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  double precision,intent(in)::rhs(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  integer,intent(in):: nits
  !     Local vars
  integer:: m,ii,jj,j,kk,ip,im,jm,jp,mm,rbi,offi,ig,jg,kg
  double precision:: ti,tj,tk,dxi2,dyi2,numer,deno
  double precision norm
#ifndef TWO_D
  double precision:: dzi2
#endif
#ifndef TWO_D

  flops = 13*g%imax*g%jmax*g%kmax/2 ! R/B
  ig = getIglobalx(g)
  jg = getIglobaly(g)
  kg = getIglobalz(g)

#endif
  dxi2=1.d0/g%dxg**2
  dyi2=1.d0/g%dyg**2
#ifndef TWO_D
  dzi2=1.d0/g%dzg**2
#endif
  do m=1,nits
     ! red/black
     do rbi = 0,1
#ifdef HAVE_PETSC
        call PetscLogEventBegin(events(2),ierr)
        call PetscLogFlops(flops,ierr)
#endif 
        do kk=1,g%kmax
           do jj=1,g%jmax
              offi = mod(ig+jj+jg+kk+kg-1+rbi,2)+1
              do ii=offi,g%imax,2
                 ti = dxi2*(phi(ii+1,jj,kk,1)+phi(ii-1,jj,kk,1))
                 tj = dyi2*(phi(ii,jj+1,kk,1)+phi(ii,jj-1,kk,1))
                 !     set
                 numer = rhs(ii,jj,kk,1) - ti - tj
                 deno = - 2.d0*(dxi2 + dyi2)
#ifndef TWO_D
                 ! Z direction
                 !     set
                 tk = dzi2*(phi(ii,jj,kk+1,1)+phi(ii,jj,kk-1,1))
                 numer = numer - tk
                 deno = deno - 2.d0*dzi2
#endif
                 phi(ii,jj,kk,1) = numer/deno
              enddo       ! ii
           enddo          ! jj
        enddo             ! kk
#ifdef HAVE_PETSC
        call PetscLogEventEnd(events(2),ierr)
#endif
        call SetBCs(phi,g)
     enddo                ! r/b i
  enddo                   ! iters
end subroutine GSRB_const_Lap
!-----------------------------------------------------------------------
subroutine Jacobi_const_Lap(phi,rhs,g,nits)
  use GridModule
  use mpistuff
  implicit none
  type(pe_patch):: g
  double precision:: phi(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  double precision,intent(in)::rhs(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  integer,intent(in):: nits
  !     Local vars
  double precision::Res(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  double precision::Dk(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  double precision::Aux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  double precision:: diagi,dxl,dyl,ti,tj,tk,dxi2,dyi2,dzi2,dzl
  double precision:: over,under,rhok,rhokp1,beta,alpha,delta,theta,s1,ct1,ct2,omega
  integer::k
#ifndef TWO_D
#ifdef HAVE_PETSC
  flops = 13*g%imax*g%jmax*g%kmax
  call PetscLogEventBegin(events(4),ierr)
  call PetscLogFlops(flops,ierr)
#endif
#endif
  dxl=g%dxg
  dyl=g%dyg
  dxi2=1.d0/dxl**2
  dyi2=1.d0/dyl**2
#ifndef TWO_D
  dzl=g%dzg
  dzi2=1.d0/dzl**2
#endif  
  diagi = -2.d0*(dxi2+dyi2)
#ifndef TWO_D
  diagi = -2.d0*(dxi2+dyi2+dzi2)
#endif
  diagi = 1.d0/diagi
  rhok = 2.d0 ! max eigen of D^-1A

  if (.true.) then
     omega = diagi*4.d0/(3.d0*rhok)
     do k=0,nits-1
        call Apply_const_Lap(Aux,phi,g)
        phi = phi + omega*(rhs - Aux)
     enddo 
  else
     over = 1.0d0
     under = 1.d0/2.d0     
     beta =  rhok * over 
     alpha = rhok * under
     delta = (beta - alpha)/2.
     theta = (beta + alpha)/2.
     s1 = theta / delta
     rhok = 1./s1

     call Apply_const_Lap(Dk,phi,g)
     Dk = rhs - Dk
     Dk = diagi * Dk
     ct1 = 1.d0/theta

     Dk = Dk*ct1
     phi = phi + Dk
     
     do k=0,nits-2
        rhokp1 = 1.d0/(2.d0*s1 - rhok)
        ct1 = rhokp1*rhok
        ct2 = 2.d0*rhokp1/delta
        rhok = rhokp1

        call Apply_const_Lap(Aux,phi,g)
        Aux = rhs - Aux
        Res = diagi * Aux
        
        ! Dk[] = rhokp1 * rhok * Dk[] + 2. * rhokp1 * res[] / ( delta * diag[] )
        Dk =  ct1*Dk + ct2*Res
        phi = phi + Dk
     enddo                   ! iters
  end if
end subroutine Jacobi_const_Lap

!-----------------------------------------------------------------------
subroutine Apply_const_Lap(uxo,ux,g)
  use GridModule
  use mpistuff
  use pms, only:verbose
  implicit none
  type(pe_patch):: g
  double precision,intent(out):: &
       uxo(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  double precision,intent(in):: &
       ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,1)
  !     Local vars
  integer its
  integer:: ii,jj,kk
  double precision:: dxl,dyl,ti,tj,tk,dxi2,dyi2,dzi2,dzl
  double precision norm
#ifndef TWO_D
#ifdef HAVE_PETSC
  flops = 13*g%imax*g%jmax*g%kmax
  call PetscLogEventBegin(events(4),ierr)
  call PetscLogFlops(flops,ierr)
#endif
#endif
  dxl=g%dxg
  dyl=g%dyg
  dxi2=1.d0/dxl**2
  dyi2=1.d0/dyl**2
#ifndef TWO_D
  dzl=g%dzg
  dzi2=1.d0/dzl**2
#endif  
  do kk=1,g%kmax,1
     do jj=1,g%jmax,1
        do ii=1,g%imax,1                    
           ti =   dxi2*(ux(ii+1,jj,kk,1)+ux(ii-1,jj,kk,1))&
                + dyi2*(ux(ii,jj+1,kk,1)+ux(ii,jj-1,kk,1))&
                - 2.d0*(dxi2+dyi2)*ux(ii,jj,kk,1)
#ifndef TWO_D
           ti = ti + dzi2*(ux(ii,jj,kk+1,1)+ux(ii,jj,kk-1,1))&
                - 2.d0*dzi2*ux(ii,jj,kk,1)
#endif
           uxo(ii,jj,kk,1) = ti               
        enddo               ! ii
     enddo                  ! jj
  enddo                     ! kk
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(4),ierr)
#endif
  call SetBCs(uxo,g)

  if (verbose.gt.2) then
     tk = norm(uxo,g,2); if(mype==0)write(6,'(A,E14.6)')'        Apply_const_Lap: done |u|=',tk
  end if
  return
end subroutine Apply_const_Lap
!-----------------------------------------------------------------------
subroutine formExactU(exact,g)
  use GridModule
  use mpistuff
  use pms
  use domain
  implicit none
  type(pe_patch),intent(in):: g
  double precision,intent(out)::exact(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar) 
  !
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(7),ierr)
#endif
  ig = getIglobalx(g)
  jg = getIglobaly(g)
  kg = getIglobalz(g)
  ! have to do exact by hand - valid region plus ghost
  do kk=g%kvallo-nsg,g%kvalhi+nsg
     do jj=g%jvallo-nsg,g%jvalhi+nsg
        do ii=g%ivallo-nsg,g%ivalhi+nsg
           coord(1) = xl+(ig+ii-1)*g%dxg-0.5*g%dxg
           coord(2) = yl+(jg+jj-1)*g%dyg-0.5*g%dyg
#ifndef TWO_D
           coord(3) = zl+(kg+kk-1)*g%dzg-0.5*g%dzg
#endif
           select case (problemType)
           case(0)              
              !write (6,'(A18,E14.6,E14.6,E14.6)') "formExactU: coord=",coord
              x2 = coord*coord
              a = x2*(x2-1.);              
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
end subroutine formExactU
!-----------------------------------------------------------------------
subroutine formGradError(exactu,ux,g,gerr,order)
  use GridModule
  use pms
  use mpistuff
  implicit none
  type(pe_patch),intent(in):: g
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
  do kk=g%kvallo,g%kvalhi
     do jj=g%jvallo,g%jvalhi
        do ii=g%ivallo,g%ivalhi
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
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(7),ierr)
#endif
end subroutine formGradError
!-----------------------------------------------------------------------
subroutine FormRHS(rhs,g)
  use GridModule
  use pms
  use domain
  use mpistuff
  implicit none
  type(pe_patch),intent(in):: g
  double precision,intent(out)::rhs(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  
  integer:: ii,jj,kk
  integer:: ig,jg,kg
  double precision::coord(3),x2(3),a(3),b(3)
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(7),ierr)
#endif 
  ig = getIglobalx(g)
  jg = getIglobaly(g)
  kg = getIglobalz(g)
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
