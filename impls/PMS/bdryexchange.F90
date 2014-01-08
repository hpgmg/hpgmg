!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
subroutine SetBCs(ux,g)
  use GridModule
  use mpistuff
  implicit none

  type(proc_patch),intent(in) :: g
  double precision:: ux(g%ilo:g%ihi,g%jlo:g%jhi,&
       g%klo:g%khi,nvar)
#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(6),ierr)
#endif
  call SetXBC(ux,g)
  call XExchange(ux,g)
  call SetYBC(ux,g)
  call YExchange(ux,g)
#ifndef TWO_D
  call SetZBC(ux,g)
  call ZExchange(ux,g)
#endif
#ifdef HAVE_PETSC
  call PetscLogEventEnd(events(6),ierr)
#endif
  return
end subroutine SetBCs
!-----------------------------------------------------------------------
subroutine XExchange(ux,g)
  use mpistuff
  use tags
  use GridModule
  implicit none
  
  type(proc_patch),intent(in) :: g
  double precision:: ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  double precision:: xbuffer_send(ng,g%jmax,g%kmax,nvar)
  double precision:: xbuffer_recv(ng,g%jmax,g%kmax,nvar)
	
  integer:: XBUFFSIZE
  integer:: jj,mm,kk,idest
  integer:: msg_id_send_x_low
  integer:: msg_id_send_x_hi
  integer:: msg_id_recv_x_low
  integer:: msg_id_recv_x_hi  
  !       -------X DIRECTION COMMUNICATION
  !       Update x-low boundaries
  XBUFFSIZE = nvar*(ng)*(g%jmax)*(g%kmax)

#ifndef XPERIODIC
  if (g%iprocx .gt. 1) then
#endif
     call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%left, MSG_XCH_XLOW_TAG, g%comm3D,msg_id_recv_x_low,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (g%iprocx .lt. g%nprocx) then
#endif
     do kk=1,g%kmax
        do jj = 1,g%jmax
           do mm = 1, ng
              xbuffer_send(mm,jj,kk,:) = ux(g%imax+1-mm,jj,kk,:)
           enddo
        enddo
     enddo
     if(g%right<0)stop 'g%right<0'
     call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%right, MSG_XCH_XLOW_TAG, g%comm3D,msg_id_send_x_low,&
          ierr)
     Call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (g%iprocx .gt. 1) then
#endif
     call MPI_Wait(msg_id_recv_x_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
     do kk=1,g%kmax
        do jj = 1,g%jmax
           do mm = 1, ng
              ux(1-mm,jj,kk,:) = xbuffer_recv(mm,jj,kk,:)
           enddo
        enddo
     enddo
#ifndef XPERIODIC
  endif
#endif 
 
#ifndef XPERIODIC
  if (g%iprocx .lt. g%nprocx) then
#endif
     call MPI_Wait(msg_id_send_x_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
  endif
#endif  
  !	update x-high boundaries
#ifndef XPERIODIC
  if (g%iprocx .lt. g%nprocx) then
#endif
     call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%right, MSG_XCH_XHI_TAG, g%comm3D,msg_id_recv_x_hi,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (g%iprocx .gt. 1) then
#endif
     do kk=1,g%kmax
        do jj = 1,g%jmax
           do mm = 1, ng
              xbuffer_send(mm,jj,kk,:) = ux(mm,jj,kk,:)
           enddo
        enddo
     enddo
     if(g%left<0)stop 'g%left<0'
     call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%left, MSG_XCH_XHI_TAG, g%comm3D, msg_id_send_x_hi,ierr)
     Call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (g%iprocx .lt. g%nprocx) then
#endif
     call MPI_Wait(msg_id_recv_x_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
     do kk=1,g%kmax
        do jj = 1,g%jmax
           do mm = 1, ng
              ux(g%imax+mm,jj,kk,:) = xbuffer_recv(mm,jj,kk,:)
           enddo
        enddo
     enddo
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (g%iprocx .gt. 1) then
#endif
     call MPI_Wait(msg_id_send_x_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
  endif
#endif
  return
end subroutine XExchange

!-----------------------------------------------------------------------
subroutine YExchange(ux,g)
  use mpistuff
  use tags
  use GridModule
  implicit none

  type(proc_patch),intent(in) :: g
  double precision:: ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  double precision:: ybuffer_send(g%ilo:g%ihi,ng,g%kmax,nvar)
  double precision:: ybuffer_recv(g%ilo:g%ihi,ng,g%kmax,nvar)
  integer:: YBUFFSIZE
  integer:: ii,mm,kk,idest
  integer:: msg_id_send_y_low
  integer:: msg_id_send_y_hi
  integer:: msg_id_recv_y_low
  integer:: msg_id_recv_y_hi
  ! -------Y DIRECTION COMMUNICATION
  !	update y-low boundaries
  YBUFFSIZE=(g%ihi-g%ilo+1)*g%kmax*nvar*ng

#ifndef YPERIODIC
  if (g%iprocy .gt. 1) then
#endif  
     call  MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%bottom,MSG_XCH_YLOW_TAG, g%comm3D,msg_id_recv_y_low,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
  endif
#endif  
  
#ifndef YPERIODIC
  if (g%iprocy .lt. g%nprocy) then
#endif  
     do kk=1,g%kmax
        do ii = g%ilo,g%ihi
           do mm = 1, ng
              ybuffer_send(ii,mm,kk,:) = ux(ii,g%jmax+1-mm,kk,:)
           enddo
        enddo
     enddo
     if(g%top<0)stop 'g%top<0'
     call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%top, MSG_XCH_YLOW_TAG, g%comm3D, msg_id_send_y_low, ierr)
     Call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
  endif
#endif

#ifndef YPERIODIC
  if(g%iprocy .gt. 1) then
#endif  
     call MPI_Wait(msg_id_recv_y_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)     
     do kk=1,g%kmax
        do ii = g%ilo,g%ihi
           do mm = 1, ng
              ux(ii,1-mm,kk,:) = ybuffer_recv(ii,mm,kk,:)
           enddo
        enddo
     enddo
#ifndef YPERIODIC
  endif
#endif  
  
#ifndef YPERIODIC
  if (g%iprocy .lt. g%nprocy) then
#endif
     call MPI_Wait(msg_id_send_y_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
  endif
#endif
  !	update y-high boundaries
#ifndef YPERIODIC
  if (g%iprocy .lt. g%nprocy) then
#endif  
     call MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%top, MSG_XCH_YHI_TAG, g%comm3D, msg_id_recv_y_hi,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
  endif
#endif  
  
#ifndef YPERIODIC
  if (g%iprocy .gt. 1) then
#endif
     do kk=1,g%kmax
        do mm = 1, ng
           do ii = g%ilo,g%ihi
              ybuffer_send(ii,mm,kk,:) = ux(ii,mm,kk,:)
           enddo
        enddo
     enddo
     if(g%bottom<0)stop 'g%bottom<0'
     call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%bottom, MSG_XCH_YHI_TAG, g%comm3D, msg_id_send_y_hi,ierr)  
     call errorhandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
  endif
#endif

#ifndef YPERIODIC
  if (g%iprocy .lt. g%nprocy) then
#endif
     call MPI_Wait(msg_id_recv_y_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
     do kk=1,g%kmax
        do mm = 1, ng
           do ii = g%ilo,g%ihi
              ux(ii,g%jmax+mm,kk,:) = ybuffer_recv(ii,mm,kk,:)
           enddo
        enddo
     enddo
#ifndef YPERIODIC
  endif
#endif    
#ifndef YPERIODIC
  if (g%iprocy .gt. 1) then
#endif  
     call MPI_Wait(msg_id_send_y_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
  endif
#endif  	
  return
end subroutine YExchange

!-----------------------------------------------------------------------
subroutine ZExchange(ux,g)
  use mpistuff
  use GridModule
  use tags
  implicit none
  
  type(proc_patch),intent(in) :: g
  double precision:: ux(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  double precision:: zbuffer_send(g%ilo:g%ihi,g%jlo:g%jhi,ng,nvar)
  double precision:: zbuffer_recv(g%ilo:g%ihi,g%jlo:g%jhi,ng,nvar)
  integer:: ZBUFFSIZE
  integer:: ii,jj,mm,idest
  integer:: msg_id_send_z_low
  integer:: msg_id_send_z_hi
  integer:: msg_id_recv_z_low
  integer:: msg_id_recv_z_hi
  !       -------Z DIRECTION COMMUNICATION
  !	update z-low boundaries
  ZBUFFSIZE=(g%ihi-g%ilo+1)*(g%jhi-g%jlo+1)*nvar*ng

#ifndef ZPERIODIC
  if ( g%iprocz > 1 ) then
#endif
     call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%behind, MSG_XCH_ZHI_TAG, g%comm3D, msg_id_recv_z_hi, ierr) 
     call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
  endif
#endif  

#ifndef ZPERIODIC
  if ( g%iprocz < g%nprocz ) then
#endif  
     do mm = 1,ng
        do jj = g%jlo, g%jhi
           do ii = g%ilo,g%ihi
              zbuffer_send(ii,jj,mm,:) = ux(ii,jj,g%kmax+1-mm,:)
           enddo
        enddo
     enddo
     if(g%forward<0)stop 'g%forward<0'
     call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%forward, MSG_XCH_ZHI_TAG, g%comm3D, msg_id_send_z_hi, ierr)
     call ErrorHandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
  endif
#endif  
  
#ifndef ZPERIODIC
  if ( g%iprocz > 1 ) then
#endif  
     call MPI_Wait(msg_id_recv_z_hi, status, ierr)
     call ErrorHandler(ierr,ERROR_WAIT)
     do mm = 1,ng
        do jj = g%jlo, g%jhi
           do ii = g%ilo,g%ihi
              ux(ii,jj,1-mm,:) = zbuffer_recv(ii,jj,mm,:)
           enddo
        enddo
     enddo
#ifndef ZPERIODIC
  endif
#endif  

#ifndef ZPERIODIC
  if ( g%iprocz < g%nprocz ) then
#endif  
     call MPI_Wait(msg_id_send_z_hi, status, ierr)
     call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
  endif
#endif       
  !       update z-high boundaries
#ifndef ZPERIODIC
  if ( g%iprocz < g%nprocz ) then
#endif  
     call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%forward, MSG_XCH_ZLOW_TAG, g%comm3D,msg_id_recv_z_low, ierr)
     call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
  endif
#endif  

#ifndef ZPERIODIC
  if ( g%iprocz > 1 ) then
#endif  
     do mm = 1,ng
        do jj = g%jlo, g%jhi
           do ii = g%ilo,g%ihi
              zbuffer_send(ii,jj,mm,:) = ux(ii,jj,mm,:)
           enddo
        enddo
     enddo
     if(g%behind<0)stop 'g%behind<0'
     call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          g%behind, MSG_XCH_ZLOW_TAG, g%comm3D, msg_id_send_z_low, ierr)
     call ErrorHandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
  endif
#endif 
 
#ifndef ZPERIODIC
  if ( g%iprocz < g%nprocz ) then
#endif  
     call MPI_Wait(msg_id_recv_z_low, status, ierr)
     call ErrorHandler(ierr,ERROR_WAIT)
     do mm = 1,ng
        do jj = g%jlo, g%jhi
           do ii = g%ilo,g%ihi
              ux(ii,jj,g%kmax+mm,:) = zbuffer_recv(ii,jj,mm,:)
           enddo
        enddo
     enddo
#ifndef ZPERIODIC
  endif
#endif  

#ifndef ZPERIODIC
  if ( g%iprocz > 1 ) then
#endif  
     call MPI_Wait(msg_id_send_z_low, status, ierr)
     call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
  endif
#endif  
  return
end subroutine ZExchange

!-----------------------------------------------------------------------
subroutine ErrorHandler(mpierr,errortype)
  use mpistuff
  implicit none
  integer::mpierr,errortype
  if(mpierr.ne.MPI_SUCCESS) then
     write(0,*) 'PMS: MPI RETURN VALUE',mype,mpierr,errortype
  endif
  return
end subroutine ErrorHandler
!-----------------------------------------------------------------
subroutine SetXBC(u,g)
  use GridModule
  use domain
  implicit none
  type(proc_patch),intent(in):: g
  double precision::u(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  
  integer xbdry_type
  integer:: ii,jj,mm,kk
  
  xbdry_type = 1 ! 0: neumann; 1: diri -- could switch on type
  
#ifndef XPERIODIC
  if (g%iprocx .eq. 1) then
     ii=0
     if(xbdry_type.eq.0) then !zero grad
        do kk=1,g%kmax,1
           do jj=1,g%jmax,1
              do mm=0,ng-1,1
                 u(ii-mm,jj,kk,:)=u(ii+mm+1,jj,kk,:)
              enddo
           enddo
        enddo
     else if(xbdry_type.eq.1) then ! diri
        do kk=1,g%kmax,1
           do jj=1,g%jmax,1
              do mm=0,ng-1,1
                 u(ii-mm,jj,kk,1)=-u(ii+mm+1,jj,kk,1)
              enddo
           enddo
        enddo
     endif
  endif  
  !  If bdry_type then typeing boundary
  if (g%iprocx .eq. g%nprocx) then
     ii=g%imax+1
     if(xbdry_type.eq.0) then ! zero gradient
        do kk=1,g%kmax,1
           do jj=1,g%jmax,1
              do mm=0,ng-1,1
                 u(ii+mm,jj,kk,:)=u(ii-mm-1,jj,kk,:)
              enddo
           enddo
        enddo
     else if(xbdry_type.eq.1) then ! Reflecting - perfect conductor      
        do kk=1,g%kmax,1
           do jj=1,g%jmax,1
              do mm=0,ng-1,1
                 u(ii+mm,jj,kk,1)=-u(ii-mm-1,jj,kk,1)
              enddo
           enddo
        enddo
     endif
  endif
#endif  
  return
end subroutine SetXBC
!-----------------------------------------------------------------
subroutine SetYBC(u,g)
  use GridModule
  use domain
  implicit none
  type(proc_patch),intent(in):: g
  double precision::u(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
  
  integer ybdry_type
  integer:: ii,jj,mm,kk
  
  ybdry_type = 1
  !	yl Boundary: Typeing
#ifndef YPERIODIC	
  do mm=0,ng-1,1
     if (g%iprocy .eq. 1) then
        jj=0
        if(ybdry_type.eq.0) then
           do kk=1,g%kmax
              do ii=g%ilo,g%ihi,1
                 u(ii,jj-mm,kk,:)=u(ii,jj+1+mm,kk,:)
              enddo
           enddo
        else
           do kk=1,g%kmax
              do ii=g%ilo,g%ihi,1
                 u(ii,jj-mm,kk,1)=-u(ii,jj+mm+1,kk,1)
              enddo
           enddo
        endif
     endif
     !	yr Boundary: Typeing
     if (g%iprocy .eq. g%nprocy) then
        jj=g%jmax+1
        if(ybdry_type.eq.0) then
           do kk=1,g%kmax
              do ii=g%ilo,g%ihi,1
                 u(ii,jj+mm,kk,:)=u(ii,jj-mm-1,kk,:)
              enddo
           enddo
        else
           do kk=1,g%kmax
              do ii=g%ilo,g%ihi,1
                 u(ii,jj+mm,kk,1)=-u(ii,jj-mm-1,kk,1)
              enddo
           enddo
        endif
     endif
  enddo
#endif  
  return
end subroutine SetYBC
!-----------------------------------------------------------------
subroutine SetZBC(u,g)
  use GridModule
  use domain
  implicit none
  type(proc_patch),intent(in):: g
  double precision::u(g%ilo:g%ihi,g%jlo:g%jhi,g%klo:g%khi,nvar)
 
  integer:: ii,jj,mm,kk
  integer zbdry_type
  
  zbdry_type = 1
  
#ifndef ZPERIODIC
  do mm=0,ng-1,1
     if (g%iprocz .eq. 1) then
        kk=0
        if(zbdry_type.eq.0) then
           do jj=g%jlo,g%jhi,1
              do ii=g%ilo,g%ihi,1
                 u(ii,jj,kk-mm,:)=u(ii,jj,kk+1+mm,:)
              enddo
           enddo
        else
           do jj=g%jlo,g%jhi,1
              do ii=g%ilo,g%ihi,1
                 u(ii,jj,kk-mm,1)=-u(ii,jj,kk+mm+1,1)
              enddo
           enddo
        endif
     endif
     !	zr Boundary: Typeing
     if (g%iprocz .eq. g%nprocz) then
        kk=g%kmax+1
        if(zbdry_type.eq.0) then
           do jj=g%jlo,g%jhi,1
              do ii=g%ilo,g%ihi,1
                 u(ii,jj,kk+mm,:)=u(ii,jj,kk-mm-1,:)
              enddo
           enddo
        else
           do jj=g%jlo,g%jhi,1
              do ii=g%ilo,g%ihi,1
                 u(ii,jj,kk+mm,1)=-u(ii,jj,kk-mm-1,1)
              enddo
           enddo
        endif
     endif
  enddo
#endif    
  return
end subroutine SetZBC
