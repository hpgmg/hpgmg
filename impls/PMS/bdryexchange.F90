!-----------------------------------------------------------------
!       Ravi Samtaney & Mark Adams
!       Copyright 2014
!-----------------------------------------------------------------
subroutine mysleep(dt)
!===============================================================================
implicit none
integer,dimension(8) :: t ! arguments for date_and_time
integer :: s1,s2,ms1,ms2  ! start and end times [ms]
real :: dt                ! desired sleep interval [ms]
!===============================================================================
! Get start time:
call date_and_time(values=t)
ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)

do ! check time:
   call date_and_time(values=t)
   ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
   if(ms2-ms1>=dt)exit
enddo
!===============================================================================
end subroutine mysleep

subroutine SetBCs(ux,p,t,ng,nv_optional)
  use discretization,only:nvar
  use base_data_module
  use mpistuff
  use bc_module,only:bc_valid
  implicit none
  integer,intent(in),optional::nv_optional
  type(ipoint),intent(in)::ng ! processor ghosts. BC, or stencil, ghost in 'nsg'
  type(topot),intent(in)::t
  type(patcht),intent(in)::p
  double precision:: ux(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,*)
  !
  integer::ii,jj,kk,nv
  
  if (present(nv_optional)) then
     nv = nv_optional
  else
     nv = nvar
  end if

#ifdef HAVE_PETSC
  call PetscLogEventBegin(events(6),ierr)
#endif
  if (t%comm3d==mpi_comm_null) then
     call SetXBC(ux,p,t,nv)
     call SetYBC(ux,p,t,nv)
#ifndef TWO_D
     call SetZBC(ux,p,t,nv)
#endif
  else
     call SetXBC(ux,p,t,nv)
     call XExchange(ux,p,t,ng%i,nv)
     call SetYBC(ux,p,t,nv)
     call YExchange(ux,p,t,ng%j,nv)
#ifndef TWO_D
     call SetZBC(ux,p,t,nv)
     call ZExchange(ux,p,t,ng%k,nv)
#endif
  end if
#ifdef HAVE_PETSC
     call PetscLogEventEnd(events(6),ierr)
#endif

  if (mype==-10) then
     ii=bc_valid%i; jj=bc_valid%j; kk=bc_valid%k
     print *,'SetBCs: AFTER kk-1 plane"',p%all%lo%i,p%all%hi%i,p%all%lo%j,p%all%hi%j,p%all%lo%k,p%all%hi%k
     print *,'      ',ux(ii-1,jj-1,kk-1,1),ux(ii,jj-1,kk-1,1),ux(ii+1,jj-1,kk-1,1)
     print *,'      ',ux(ii-1,jj,  kk-1,1),ux(ii,jj,  kk-1,1),ux(ii+1,jj,  kk-1,1)
     print *,'      ',ux(ii-1,jj+1,kk-1,1),ux(ii,jj+1,kk-1,1),ux(ii+1,jj+1,kk-1,1)
     print *,'      ','kk:'
     print *,'      ',ux(ii-1,jj-1,kk,  1),ux(ii,jj-1,kk,  1),ux(ii+1,jj-1,kk,  1)
     print *,'      ',ux(ii-1,jj,  kk,  1),ux(ii,jj,  kk,  1),ux(ii+1,jj,  kk,  1)
     print *,'      ',ux(ii-1,jj+1,kk,  1),ux(ii,jj+1,kk,  1),ux(ii+1,jj+1,kk,  1)
     print *,'      ','kk+1:'
     print *,'      ',ux(ii-1,jj-1,kk+1,1),ux(ii,jj-1,kk+1,1),ux(ii+1,jj-1,kk+1,1)
     print *,'      ',ux(ii-1,jj,  kk+1,1),ux(ii,jj,  kk+1,1),ux(ii+1,jj,  kk+1,1)
     print *,'      ',ux(ii-1,jj+1,kk+1,1),ux(ii,jj+1,kk+1,1),ux(ii+1,jj+1,kk+1,1)
  end if
  call flush(6)
  !call mysleep(80.)

  return
end subroutine SetBCs
!-----------------------------------------------------------------------
subroutine XExchange(ux,p,t,ng,nv)
  use base_data_module
  use mpistuff
  use tags
  implicit none
  integer,intent(in)::nv
  integer,intent(in)::ng
  type(patcht),intent(in)::p
  type(topot),intent(in)::t
  double precision:: ux(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nv)
  !
  double precision:: xbuffer_send(ng,p%max%hi%j,p%max%hi%k,nv)
  double precision:: xbuffer_recv(ng,p%max%hi%j,p%max%hi%k,nv)
  integer:: XBUFFSIZE
  integer:: jj,mm,kk
  integer:: msg_id_send_x_low
  integer:: msg_id_send_x_hi
  integer:: msg_id_recv_x_low
  integer:: msg_id_recv_x_hi  
  !       -------X DIRECTION COMMUNICATION
  !       Update x-low boundaries
  XBUFFSIZE = size(xbuffer_send)
  !call mpi_barrier(mpi_comm_world,ierr)
#ifndef XPERIODIC
  if (t%ipe%i.gt.1) then
#endif
     call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%left,MSG_XCH_XLOW_TAG,t%comm3D,msg_id_recv_x_low,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
  endif
#endif

#ifndef XPERIODIC
  if (t%ipe%i .lt. t%npe%i) then
#endif
     do kk=1,p%max%hi%k
        do jj=1,p%max%hi%j
           do mm=1,ng
              xbuffer_send(mm,jj,kk,:) = ux(p%max%hi%i+1-mm,jj,kk,:)
           enddo
        enddo
     enddo
     if (t%right<0)stop 'p%t%right<0'
     call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%right,MSG_XCH_XLOW_TAG,t%comm3D,msg_id_send_x_low,ierr)
     Call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (t%ipe%i .gt. 1) then
#endif
     call MPI_Wait(msg_id_recv_x_low,status,ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
     do kk=1,p%max%hi%k
        do jj=1,p%max%hi%j
           do mm=1,ng
              ux(1-mm,jj,kk,:) = xbuffer_recv(mm,jj,kk,:)
           enddo
        enddo
     enddo
#ifndef XPERIODIC
  endif
#endif 
 
#ifndef XPERIODIC
  if (t%ipe%i .lt. t%npe%i) then
#endif
     call MPI_Wait(msg_id_send_x_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
  endif
#endif  
  !	update x-high boundaries
#ifndef XPERIODIC
  if (t%ipe%i .lt. t%npe%i) then
#endif
     call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%right, MSG_XCH_XHI_TAG, t%comm3D,msg_id_recv_x_hi,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (t%ipe%i .gt. 1) then
#endif
     do kk=1,p%max%hi%k
        do jj=1,p%max%hi%j
           do mm=1,ng
              xbuffer_send(mm,jj,kk,:) = ux(mm,jj,kk,:)
           enddo
        enddo
     enddo
     if (t%left<0)stop 't%left<0'
     call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%left, MSG_XCH_XHI_TAG, t%comm3D, msg_id_send_x_hi,ierr)
     Call ErrorHandler(ierr,ERROR_SEND)
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (t%ipe%i .lt. t%npe%i) then
#endif
     call MPI_Wait(msg_id_recv_x_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
     do kk=1,p%max%hi%k
        do jj=1,p%max%hi%j
           do mm=1,ng
              ux(p%max%hi%i+mm,jj,kk,:) = xbuffer_recv(mm,jj,kk,:)
           enddo
        enddo
     enddo
#ifndef XPERIODIC
  endif
#endif
  
#ifndef XPERIODIC
  if (t%ipe%i .gt. 1) then
#endif
     call MPI_Wait(msg_id_send_x_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef XPERIODIC
  endif
#endif
  ! keep from getting mixed up 
  MSG_XCH_XLOW_TAG = MSG_XCH_XLOW_TAG+1
  MSG_XCH_XHI_TAG = MSG_XCH_XHI_TAG+1
  return
end subroutine XExchange
!-----------------------------------------------------------------------
subroutine YExchange(ux,p,t,ng,nv)
  use base_data_module
  use mpistuff
  use tags
  implicit none
  integer,intent(in)::nv
  integer,intent(in)::ng
  type(patcht),intent(in) :: p
  type(topot),intent(in)::t
  double precision:: ux(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nv)

  double precision:: ybuffer_send(p%all%lo%i:p%all%hi%i,ng,p%max%hi%k,nv)
  double precision:: ybuffer_recv(p%all%lo%i:p%all%hi%i,ng,p%max%hi%k,nv)
  integer:: YBUFFSIZE
  integer:: ii,mm,kk
  integer:: msg_id_send_y_low
  integer:: msg_id_send_y_hi
  integer:: msg_id_recv_y_low
  integer:: msg_id_recv_y_hi
  ! -------Y DIRECTION COMMUNICATION
  !	update y-low boundaries
  YBUFFSIZE=size(ybuffer_recv)

#ifndef YPERIODIC
  if (t%ipe%j .gt. 1) then
#endif  
     call  MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%bottom,MSG_XCH_YLOW_TAG, t%comm3D,msg_id_recv_y_low,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
  endif
#endif  

#ifndef YPERIODIC
  if (t%ipe%j .lt. t%npe%j) then
#endif  
     do ii=p%all%lo%i,p%all%hi%i
        do kk=1,p%max%hi%k
           do mm=1,ng
              ybuffer_send(ii,mm,kk,:) = ux(ii,p%max%hi%j+1-mm,kk,:)
           enddo
        enddo
     enddo
     if (t%top<0)stop 't%top<0'
     call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%top,MSG_XCH_YLOW_TAG,t%comm3D,msg_id_send_y_low,ierr)
     Call ErrorHandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
  endif
#endif

#ifndef YPERIODIC
  if (t%ipe%j .gt. 1) then
#endif  
     call MPI_Wait(msg_id_recv_y_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)     
     do kk=1,p%max%hi%k
        do ii=p%all%lo%i,p%all%hi%i
           do mm=1,ng
              ux(ii,1-mm,kk,:) = ybuffer_recv(ii,mm,kk,:)
           enddo
        enddo
     enddo
#ifndef YPERIODIC
  endif
#endif  
  
#ifndef YPERIODIC
  if (t%ipe%j .lt. t%npe%j) then
#endif
     call MPI_Wait(msg_id_send_y_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
  endif
#endif
  !	update y-high boundaries
#ifndef YPERIODIC
  if (t%ipe%j .lt. t%npe%j) then
#endif  
     call MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%top, MSG_XCH_YHI_TAG, t%comm3D, msg_id_recv_y_hi,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef YPERIODIC
  endif
#endif  
  
#ifndef YPERIODIC
  if (t%ipe%j .gt. 1) then
#endif
     do kk=1,p%max%hi%k
        do ii=p%all%lo%i,p%all%hi%i
           do mm=1,ng
              ybuffer_send(ii,mm,kk,:) = ux(ii,mm,kk,:)
           enddo
        enddo
     enddo
     if (t%bottom<0)stop 't%bottom<0'
     call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%bottom, MSG_XCH_YHI_TAG, t%comm3D, msg_id_send_y_hi,ierr)  
     call errorhandler(ierr,ERROR_SEND)
#ifndef YPERIODIC
  endif
#endif

#ifndef YPERIODIC
  if (t%ipe%j .lt. t%npe%j) then
#endif
     call MPI_Wait(msg_id_recv_y_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
     do kk=1,p%max%hi%k
        do ii=p%all%lo%i,p%all%hi%i
           do mm=1,ng
              ux(ii,p%max%hi%j+mm,kk,:) = ybuffer_recv(ii,mm,kk,:)
           enddo
        enddo
     enddo
#ifndef YPERIODIC
  endif
#endif    
#ifndef YPERIODIC
  if (t%ipe%j .gt. 1) then
#endif  
     call MPI_Wait(msg_id_send_y_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef YPERIODIC
  endif
#endif  
  ! keep from getting mixed up 
  MSG_XCH_YLOW_TAG = MSG_XCH_YLOW_TAG+1
  MSG_XCH_YHI_TAG = MSG_XCH_YHI_TAG+1
  return
end subroutine YExchange
!-----------------------------------------------------------------------
subroutine ZExchange(ux,p,t,ng,nv)
  use base_data_module
  use mpistuff
  use tags
  implicit none
  integer,intent(in)::nv
  integer,intent(in)::ng
  type(patcht),intent(in) :: p
  type(topot),intent(in)::t
  double precision:: ux(&
       p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nv)

  double precision:: zbuffer_send(p%all%lo%i:p%all%hi%i,p%all%lo%j:p%all%hi%j,ng,nv)
  double precision:: zbuffer_recv(p%all%lo%i:p%all%hi%i,p%all%lo%j:p%all%hi%j,ng,nv)
  integer:: ZBUFFSIZE
  integer:: ii,jj,mm
  integer:: msg_id_send_z_low
  integer:: msg_id_send_z_hi
  integer:: msg_id_recv_z_low
  integer:: msg_id_recv_z_hi
  ! -------Z DIRECTION COMMUNICATION
  !	update z-low boundaries
  ZBUFFSIZE=size(zbuffer_recv)

#ifndef ZPERIODIC
  if (t%ipe%k .gt. 1) then
#endif  
     call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%behind,MSG_XCH_ZLOW_TAG, t%comm3D,msg_id_recv_z_low,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
  endif
#endif  
  
#ifndef ZPERIODIC
  if (t%ipe%k .lt. t%npe%k) then
#endif  
     do ii=p%all%lo%i,p%all%hi%i
        do jj=p%all%lo%j,p%all%hi%j
           do mm=1,ng
              zbuffer_send(ii,jj,mm,:) = ux(ii,jj,p%max%hi%k+1-mm,:)
           enddo
        enddo
     enddo
     if (t%forward<0)stop 't%forward<0'
     call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%forward,MSG_XCH_ZLOW_TAG,t%comm3D,msg_id_send_z_low,ierr)
     Call ErrorHandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
  endif
#endif

#ifndef ZPERIODIC
  if (t%ipe%k .gt. 1) then
#endif  
     call MPI_Wait(msg_id_recv_z_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)     
     do ii=p%all%lo%i,p%all%hi%i
        do jj=p%all%lo%j,p%all%hi%j
           do mm=1,ng
              ux(ii,jj,1-mm,:) = zbuffer_recv(ii,jj,mm,:)
           enddo
        enddo
     enddo
#ifndef ZPERIODIC
  endif
#endif  
  
#ifndef ZPERIODIC
  if (t%ipe%k .lt. t%npe%k) then
#endif
     call MPI_Wait(msg_id_send_z_low, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
  endif
#endif
  !	update z-high boundaries
#ifndef ZPERIODIC
  if (t%ipe%k .lt. t%npe%k) then
#endif  
     call MPI_Irecv(zbuffer_recv, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%forward, MSG_XCH_ZHI_TAG, t%comm3D, msg_id_recv_z_hi,ierr)
     Call ErrorHandler(ierr,ERROR_RECV)
#ifndef ZPERIODIC
  endif
#endif  
  
#ifndef ZPERIODIC
  if (t%ipe%k .gt. 1) then
#endif
     do ii=p%all%lo%i,p%all%hi%i
        do jj=p%all%lo%j,p%all%hi%j
           do mm=1,ng
              zbuffer_send(ii,jj,mm,:) = ux(ii,jj,mm,:)
           enddo
        enddo
     enddo
     if (t%behind<0)stop 't%behind<0'
     call MPI_Isend(zbuffer_send, ZBUFFSIZE, MPI_DOUBLE_PRECISION,&
          t%behind, MSG_XCH_ZHI_TAG, t%comm3D, msg_id_send_z_hi,ierr)  
     call errorhandler(ierr,ERROR_SEND)
#ifndef ZPERIODIC
  endif
#endif

#ifndef ZPERIODIC
  if (t%ipe%k .lt. t%npe%k) then
#endif
     call MPI_Wait(msg_id_recv_z_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
     do ii=p%all%lo%i,p%all%hi%i
        do jj=p%all%lo%j,p%all%hi%j
           do mm=1,ng
              ux(ii,jj,p%max%hi%k+mm,:) = zbuffer_recv(ii,jj,mm,:)
           enddo
        enddo
     enddo
#ifndef ZPERIODIC
  endif
#endif    
#ifndef ZPERIODIC
  if (t%ipe%k .gt. 1) then
#endif  
     call MPI_Wait(msg_id_send_z_hi, status, ierr)
     Call ErrorHandler(ierr,ERROR_WAIT)
#ifndef ZPERIODIC
  endif
#endif  
  ! keep from getting mixed up 
  MSG_XCH_ZLOW_TAG = MSG_XCH_ZLOW_TAG+1
  MSG_XCH_ZHI_TAG = MSG_XCH_ZHI_TAG+1
  return
end subroutine ZExchange
!-----------------------------------------------------------------
subroutine SetXBC(u,p,t,nv)
  use discretization,only:nsg
  use base_data_module
  use domain
  implicit none
  integer,intent(in)::nv
  type(patcht),intent(in):: p
  type(topot),intent(in)::t
  double precision::u(p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nv)
  !
  integer xbdry_type
  integer:: ii,jj,mm,kk

  xbdry_type = 1 ! 0: neumann; 1: diri -- could switch on type
#ifndef XPERIODIC
  do mm=0,nsg%i-1,1
     if (t%ipe%i .eq. 1) then
        ii=0
        if (xbdry_type.eq.0) then !zero grad
           do kk=1,p%max%hi%k,1
              do jj=1,p%max%hi%j,1                 
                 u(ii-mm,jj,kk,:)=u(ii+mm+1,jj,kk,:)
              enddo
           enddo
        else if (xbdry_type.eq.1) then ! diri
           do kk=1,p%max%hi%k,1
              do jj=1,p%max%hi%j,1
                 u(ii-mm,jj,kk,:)=-u(ii+mm+1,jj,kk,:)
              enddo
           enddo
        end if
     end if
     ! If bdry_type then typeing boundary
     if (t%ipe%i .eq. t%npe%i) then
        ii=p%max%hi%i+1
        if (xbdry_type.eq.0) then ! zero gradient
           do kk=1,p%max%hi%k,1
              do jj=1,p%max%hi%j,1
                 u(ii+mm,jj,kk,:)=u(ii-mm-1,jj,kk,:)
              enddo
           enddo
        else if (xbdry_type.eq.1) then ! Reflecting - perfect conductor      
           do kk=1,p%max%hi%k,1
              do jj=1,p%max%hi%j,1
                 u(ii+mm,jj,kk,:)=-u(ii-mm-1,jj,kk,:)
              enddo
           enddo
        end if
     endif
  end do
#endif
  return
end subroutine SetXBC
!-----------------------------------------------------------------
subroutine SetYBC(u,p,t,nv)
  use discretization,only:nsg
  use base_data_module
  use domain
  implicit none
  integer,intent(in)::nv
  type(patcht),intent(in):: p
  double precision::u(p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nv)
  type(topot),intent(in)::t

  integer ybdry_type
  integer:: ii,jj,mm,kk
  
  ybdry_type = 1
  !	yl Boundary: Typeing
#ifndef YPERIODIC	
  do mm=0,nsg%j-1,1
     if (t%ipe%j .eq. 1) then
        jj=0
        if (ybdry_type.eq.0) then
           do kk=1,p%max%hi%k
              do ii=p%all%lo%i,p%all%hi%i ! this pet corners but they are wrong!!
                 u(ii,jj-mm,kk,:)=u(ii,jj+1+mm,kk,:)
              enddo
           enddo
        else
           do kk=1,p%max%hi%k
              do ii=p%all%lo%i,p%all%hi%i
                 u(ii,jj-mm,kk,:)=-u(ii,jj+mm+1,kk,:)
              enddo
           enddo
        endif
     endif
     !	yr Boundary: Typeing
     if (t%ipe%j .eq. t%npe%j) then
        jj=p%max%hi%j+1
        if (ybdry_type.eq.0) then
           do kk=1,p%max%hi%k
              do ii=p%all%lo%i,p%all%hi%i
                 u(ii,jj+mm,kk,:)=u(ii,jj-mm-1,kk,:)
              enddo
           enddo
        else
           do kk=1,p%max%hi%k
              do ii=p%all%lo%i,p%all%hi%i
                 u(ii,jj+mm,kk,:)=-u(ii,jj-mm-1,kk,:)
              enddo
           enddo
        endif
     endif
  enddo
#endif  
  return
end subroutine SetYBC
!-----------------------------------------------------------------
subroutine SetZBC(u,p,t,nv)
  use discretization,only:nsg
  use base_data_module
  use domain
  implicit none
  integer,intent(in)::nv
  type(patcht),intent(in):: p
  double precision::u(p%all%lo%i:p%all%hi%i,&
       p%all%lo%j:p%all%hi%j,&
       p%all%lo%k:p%all%hi%k,nv)
  type(topot),intent(in)::t

  integer:: ii,jj,mm,kk
  integer zbdry_type
  
  zbdry_type = 1
#ifndef ZPERIODIC
  do mm=0,nsg%k-1,1
     if (t%ipe%k .eq. 1) then
        kk=0
        if (zbdry_type.eq.0) then
           do ii=p%all%lo%i,p%all%hi%i
              do jj=p%all%lo%j,p%all%hi%j
                 u(ii,jj,kk-mm,:)=u(ii,jj,kk+1+mm,:)
              enddo
           enddo
        else
           do ii=p%all%lo%i,p%all%hi%i
              do jj=p%all%lo%j,p%all%hi%j
                 u(ii,jj,kk-mm,:)=-u(ii,jj,kk+mm+1,:)
              enddo
           enddo
        endif
     endif
     !	zr Boundary: Typeing
     if (t%ipe%k .eq. t%npe%k) then
        kk=p%max%hi%k+1
        if (zbdry_type.eq.0) then
           do ii=p%all%lo%i,p%all%hi%i
              do jj=p%all%lo%j,p%all%hi%j
                 u(ii,jj,kk+mm,:)=u(ii,jj,kk-mm-1,:)
              enddo
           enddo
        else
           do ii=p%all%lo%i,p%all%hi%i
              do jj=p%all%lo%j,p%all%hi%j
                 u(ii,jj,kk+mm,:)=-u(ii,jj,kk-mm-1,:)
              enddo
           enddo
        endif
     endif
  enddo
#endif    
  return
end subroutine SetZBC
!-----------------------------------------------------------------
subroutine ErrorHandler(mpierr,errortype)
  use mpistuff
  implicit none
  integer::mpierr,errortype
  if (mpierr.ne.MPI_SUCCESS) then
     write(0,*) 'PMS: MPI RETURN VALUE',mype,mpierr,errortype
  endif
  return
end subroutine ErrorHandler
