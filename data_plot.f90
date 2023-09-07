subroutine data_plot

use mpi
use para
use euler
use lagrange

implicit none
integer::nsend,nrecv,id,ib,ie,jb,je,sendindex,recvindex
integer::status(MPI_STATUS_SIZE),tag,request1,request2,request3,request4
double precision:: pnorm8,pnorm
double precision,dimension(:),allocatable::xcc,ycc
double precision,dimension(:,:),allocatable::uu,vv,ww,pp
double precision,dimension(:),allocatable::usend,urecv,vsend,vrecv
double precision,dimension(:),allocatable::psend,precv,wsend,wrecv
! collect u,v,w,p,xc,yc
! interpolate uv
! store data

request1=MPI_REQUEST_NULL
request2=MPI_REQUEST_NULL
request3=MPI_REQUEST_NULL
request4=MPI_REQUEST_NULL
tag=1234

! interpolate u,v at center grid points
!call interpolateuv

! initialize u,v,p,vorticity on center grid points
if(myid.eq.0) then
  allocate(uu(1:nx,1:ny))
  allocate(vv(1:nx,1:ny))
  allocate(pp(1:nx,1:ny))
  allocate(ww(1:nx,1:ny))
endif

! initialize send buffer
nsend=(jpend-jpbeg+1)*(ipend-ipbeg+1)
allocate(usend(1:nsend))
allocate(vsend(1:nsend))
allocate(psend(1:nsend))
allocate(wsend(1:nsend))
usend=0
vsend=0
psend=0
wsend=0
do j=jpbeg,jpend
  do i=ipbeg,ipend
    sendindex=(i-ipbeg+1)+(j-jpbeg)*(ipend-ipbeg+1)
    usend(sendindex)=myid+1
!    vsend(sendindex)=vcc(i,j)
    !psend(sendindex)=p(i,j)
    !wsend(sendindex)=wo(i,j)
  enddo
enddo
!write(*,*),'send',myid,usend(1),usend(nsend)   
!write(*,*),myid,ipbeg,ipend,jpbeg,jpend,nsend

! send u,v,p,w to processor 0
call MPI_Isend(usend,nsend,MPI_DOUBLE_PRECISION,&
                  0,tag,comm2d,request1,ierr)
!call MPI_Isend(vsend,nsend,MPI_DOUBLE_PRECISION,&
!                 0,tag,comm2d,request2,ierr)
!call MPI_Isend(psend,nsend,MPI_DOUBLE_PRECISION,&
!                 0,tag,comm2d,request3,ierr)
!call MPI_Isend(wsend,nsend,MPI_DOUBLE_PRECISION,&
!                 0,tag,comm2d,request4,ierr)

! receive u,v,p,w to processor 0
if(myid.eq.0) then
!  do id=0,nprocs-1
   do id=0,1
    ib=ProcsInfo(4,id)
    ie=ProcsInfo(5,id) 
    jb=ProcsInfo(6,id)
    je=ProcsInfo(7,id)
    nrecv=(je-jb+1)*(ie-ib+1)

! write(*,*),id,ib,ie,jb,je
   ! allocate receive buffer
    allocate(urecv(1:nrecv))
!    allocate(vrecv(1:nrecv))
!    allocate(precv(1:nrecv))
!    allocate(wrecv(1:nrecv))   
    urecv=0
    ! receive
    call MPI_Irecv(urecv,nrecv,MPI_DOUBLE_PRECISION,id,&
                      tag,comm2d,request1,ierr)
!    call MPI_Irecv(vrecv,nrecv,MPI_DOUBLE_PRECISION,id,&
!                     tag,comm2d,request2,ierr)
!    call MPI_Irecv(precv,nrecv,MPI_DOUBLE_PRECISION,id,&
!                     tag,comm2d,request3,ierr)
!    call MPI_Irecv(wrecv,nrecv,MPI_DOUBLE_PRECISION,id,&
!                     tag,comm2d,request4,ierr)
write(*,*),'recv',id,urecv(1),urecv(nrecv)
!   allocate receive information
    do j=jb,je
      do i=ib,ie
        recvindex=(i-ib+1)+(j-jb)*(ie-ib+1)
        uu(i,j)=urecv(recvindex)
!        vv(i,j)=vrecv(recvindex)
!        pp(i,j)=precv(recvindex)
!        ww(i,j)=wrecv(recvindex)
      enddo
    enddo
!    write(*,*),uu(ib,jb),uu(ie,jb),uu(ib,je),uu(ie,je)
     deallocate(urecv)
!    deallocate(vrecv)
!    deallocate(precv)
!    deallocate(wrecv)
  write(*,*),'finish 0'
  enddo !end id=0,nprocs-1
endif   
! wait send recv to finish
!call MPI_Wait(request1,status,ierr)
!call MPI_Wait(request2,status,ierr)
!call MPI_Wait(request3,status,ierr)
!call MPI_Wait(request4,status,ierr)

! output 
if(myid.eq.0) then
! streamfunction
  open(unit=11,file='OUTPUT/xc.dat',status='unknown')
  open(unit=12,file='OUTPUT/yc.dat',status='unknown')
  do j=jpbeg,jpend
    do i=ipbeg,ipend
      
    enddo
  enddo

  close(11)
  close(12)

      open(unit=38,file='OUTPUT/p.dat',status='unknown')
      open(unit=48,file='OUTPUT/d.dat',status='unknown')
      open(unit=58,file='OUTPUT/uc.dat',status='unknown')
      open(unit=68,file='OUTPUT/vc.dat',status='unknown')
      open(unit=78,file='OUTPUT/wo.dat',status='unknown')
      open(unit=88,file='OUTPUT/ph.dat',status='unknown')
      do j=1,ny
!        write(38,400)(p(i,j),i=1,nx)
!        write(48,400)(d(i,j),i=1,nx)
        write(58,400)(uu(i,j),i=1,nx)
        write(68,400)(vv(i,j),i=1,nx)
!        write(78,400)(o(i,j),i=1,nx)
!        write(88,400)(-ph(i,j),i=1,nx)
      enddo
      close(38)
      close(48)
      close(58)
      close(68)
      close(78)
      close(88)

400   format(1x,2000e16.6)
endif !end myid=0
pnorm8=ZERO
 do j = jpbeg,jpend
    do i = ipbeg,ipend
      pnorm8=max(pnorm8,abs(p(i,j)-sin(xc(i))*sin(yc(j))))
    end do
  end do
call MPI_Reduce(pnorm8,pnorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)

   if ( myid .eq. 0 ) then
     print *
!     print '(A,I2)', " Iterations = ",iters 
!     print '(A,ES16.8)',"Relative Residual Norm = ",resid
!     write(*,*) "pnorm = ",pnorm
     print *
   endif


end subroutine data_plot
