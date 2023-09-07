subroutine exchangeobject
use para
use Euler
use Lagrange
use MPI
integer:: status(MPI_STATUS_SIZE)
integer:: countobj,countvertex,countpanel,countsend,countrecv
integer:: ntotalobj4send,ntotalobj4recv,ntotalvertex4recv,ntotalpanel4recv
integer:: oldnobj4proc,nobj4delete,nobj4repeat,nlocvertex,nlocpanel
integer:: myobj,myvertex,mypanel,myproc,neighbor

integer, dimension(0:26):: nobj4send,nobj4recv
integer, dimension(:), allocatable:: locindex4sendobj,nvertex4recv,npanel4recv
integer, dimension(:), allocatable:: oldobj4proc,oldnvertex4proc,oldnpanel4proc
integer, dimension(:,:), allocatable:: oldproc4obj,obj4send,obj4recv
integer, dimension(:,:), allocatable:: oldpanel,newpanel

double precision, dimension(1:3):: A,B,C,corner0,corner1
double precision, dimension(:,:), allocatable:: newvertex,newcurvature
double precision, dimension(:,:), allocatable:: oldvertex,oldcurvature

interface
  function panelinPE(A,B,C,corner0,corner1)
    double precision, dimension(1:3), intent(in):: A,B,C,corner0,corner1
    logical panelinPE
  end function panelinPE
end interface

! determine neighboring processors for each object
allocate(oldproc4obj(0:26,nobj4proc))
oldproc4obj=proc4obj
proc4obj=0
do myobj=1,nobj4proc
  do neighbor=0,26
    myproc=nearproc(neighbor)
    do mypanel=npanel4proc(myobj-1)+1,npanel4proc(myobj)
      A=vertex(:,nvertex4proc(myobj-1)+panel(1,mypanel))
      B=vertex(:,nvertex4proc(myobj-1)+panel(2,mypanel))
      C=vertex(:,nvertex4proc(myobj-1)+panel(3,mypanel))
      corner0=allcorner0(:,myproc)
      corner1=allcorner1(:,myproc)
      if(panelinPE(A,B,C,corner0,corner1)) then
        proc4obj(neighbor,myobj)=1
        exit ! regard the object on myproc if one of its panels on myproc
      endif
    end do
  end do
end do
proc4obj=proc4obj-oldproc4obj ! delte(@ 0): -1=0-1, send (@1-26): 1=1-0, no action: 0=1-1=0-0
deallocate(oldproc4obj)

! determine number of objects out of processor myid (@0)
nobj4delete=0
do myobj=1,nobj4proc
  nobj4delete=nobj4delete-proc4obj(0,myobj)
end do

! determine number of objects to be sent to each neighbor processor
nobj4send=0
do neighbor=1,26
  nobj4send(neighbor)=0
  do myobj=1,nobj4proc
    if(proc4obj(neighbor,myobj)==1) then
      nobj4send(neighbor)=nobj4send(neighbor)+1
    end if
  end do
end do

! send and receive number of objects for exchange
nobj4recv=0
do neighbor=2,26,2 !send from one direction and receive from its opposite direction, then alternate
  call MPI_SENDRECV(nobj4send(neighbor),1,MPI_INTEGER,nearproc(neighbor),0,&
                    nobj4recv(neighbor-1),1,MPI_INTEGER,nearproc(neighbor-1),0,comm3d,status,ierr)
  call MPI_SENDRECV(nobj4send(neighbor-1),1,MPI_INTEGER,nearproc(neighbor-1),0,&
                    nobj4recv(neighbor),1,MPI_INTEGER,nearproc(neighbor),0,comm3d,status,ierr)
end do

! determine ending index of each object in object stack
do neighbor=1,26
  nobj4send(neighbor)=nobj4send(neighbor-1)+nobj4send(neighbor)
end do
ntotalobj4send=nobj4send(26)
do neighbor=1,26
  nobj4recv(neighbor)=nobj4recv(neighbor-1)+nobj4recv(neighbor)
end do
ntotalobj4recv=nobj4recv(26)

! pack objects4send in stack
allocate(obj4send(3,ntotalobj4send))
obj4send=0
allocate(locindex4sendobj(ntotalobj4send))
locindex4sendobj=0

countobj=1
do neighbor=1,26
  do myobj=1,nobj4proc
    if(proc4obj(neighbor,myobj)==1) then
      locindex4sendobj(countobj)=myobj ! local object index
      obj4send(1,countobj)=obj4proc(myobj) ! global object index
      obj4send(2,countobj)=nvertex4proc(myobj)-nvertex4proc(myobj-1) ! number of vertices
      obj4send(3,countobj)=npanel4proc(myobj)-npanel4proc(myobj-1) ! number of panels
      countobj=countobj+1
    end if
  end do
end do

! exchange object info with neighbors
allocate(obj4recv(3,ntotalobj4recv))
obj4recv=0

do neighbor=2,26,2
  countsend=nobj4send(neighbor)-nobj4send(neighbor-1)
  if(countsend/=0) then
    call MPI_SEND(obj4send(1,nobj4send(neighbor-1)+1),3*countsend,MPI_INTEGER,nearproc(neighbor),0,comm3d,status,ierr)
  end if
  countrecv=nobj4recv(neighbor-1)-nobj4recv(neighbor-2)
  if(countrecv/=0) then
    call MPI_RECV(obj4recv(1,nobj4recv(neighbor-2)+1),3*countrecv,MPI_INTEGER,nearproc(neighbor-1),0,comm3d,status,ierr)
  end if
  countsend=nobj4send(neighbor-1)-nobj4send(neighbor-2)
  if(countsend/=0) then
    call MPI_SEND(obj4send(1,nobj4send(neighbor-2)+1),3*countsend,MPI_INTEGER,nearproc(neighbor-1),0,comm3d,status,ierr)
  end if
  countrecv=nobj4recv(neighbor)-nobj4recv(neighbor-1)
  if(countrecv/=0) then
    call MPI_RECV(obj4recv(1,nobj4recv(neighbor-1)+1),3*countrecv,MPI_INTEGER,nearproc(neighbor),0,comm3d,status,ierr)
  end if
end do

deallocate(obj4send)

! determine number of panels and vertices to be received
allocate(nvertex4recv(0:ntotalobj4recv))
allocate(npanel4recv(0:ntotalobj4recv))
nvertex4recv=0
npanel4recv=0

do neighbor=1,26
  do countobj=nobj4recv(neighbor-1)+1,nobj4recv(neighbor)
    nvertex4recv(countobj)=nvertex4recv(countobj-1)+obj4recv(2,countobj)
    npanel4recv(countobj)=npanel4recv(countobj-1)+obj4recv(3,countobj)
  end do
end do
ntotalvertex4recv=nvertex4recv(ntotalobj4recv)
ntotalpanel4recv=npanel4recv(ntotalobj4recv)

! exchange vertices, panels and curvatures
allocate(newvertex(3,ntotalvertex4recv))
allocate(newpanel(3,ntotalpanel4recv))
allocate(newcurvature(3,ntotalpanel4recv))
! vertices
do neighbor=2,26,2
  do countobj=nobj4send(neighbor-1)+1,nobj4send(neighbor)
    myobj=locindex4sendobj(countobj)
    myvertex=nvertex4proc(myobj-1)+1
    countvertex=nvertex4proc(myobj)-nvertex4proc(myobj-1)
    if(countvertex/=0) then
      call MPI_SEND(vertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4recv(neighbor-2)+1,nobj4recv(neighbor-1)
    myvertex=nvertex4recv(countobj-1)+1
    countvertex=nvertex4recv(countobj)-nvertex4recv(countobj-1)
    if(countvertex/=0) then
      call MPI_RECV(newvertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor-1),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4send(neighbor-2)+1,nobj4send(neighbor-1)
    myobj=locindex4sendobj(countobj)
    myvertex=nvertex4proc(myobj-1)+1
    countvertex=nvertex4proc(myobj)-nvertex4proc(myobj-1)
    if(countvertex/=0) then
      call MPI_SEND(vertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor-1),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4recv(neighbor-1)+1,nobj4recv(neighbor)
    myvertex=nvertex4recv(countobj-1)+1
    countvertex=nvertex4recv(countobj)-nvertex4recv(countobj-1)
    if(countvertex/=0) then
      call MPI_RECV(newvertex(1,myvertex),3*countvertex,MPI_DOUBLE_PRECISION,nearproc(neighbor),0,comm3d,status,ierr)
    end if
  end do
end do
! panels
do neighbor=2,26,2
  do countobj=nobj4send(neighbor-1)+1,nobj4send(neighbor)
    myobj=locindex4sendobj(countobj)
    mypanel=npanel4proc(myobj-1)+1
    countpanel=npanel4proc(myobj)-npanel4proc(myobj-1)
    if(countpanel/=0) then
      call MPI_SEND(panel(1,mypanel),3*countpanel,MPI_INTEGER,nearproc(neighbor),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4recv(neighbor-2)+1,nobj4recv(neighbor-1)
    mypanel=npanel4recv(countobj-1)+1
    countpanel=npanel4recv(countobj)-npanel4recv(countobj-1)
    if(countpanel/=0) then
      call MPI_RECV(newpanel(1,mypanel),3*countpanel,MPI_INTEGER,nearproc(neighbor-1),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4send(neighbor-2)+1,nobj4send(neighbor-1)
    myobj=locindex4sendobj(countobj)
    mypanel=npanel4proc(myobj-1)+1
    countpanel=npanel4proc(myobj)-npanel4proc(myobj-1)
    if(countpanel/=0) then
      call MPI_SEND(panel(1,mypanel),3*countpanel,MPI_INTEGER,nearproc(neighbor-1),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4recv(neighbor-1)+1,nobj4recv(neighbor)
    mypanel=npanel4recv(countobj-1)+1
    countpanel=npanel4recv(countobj)-npanel4recv(countobj-1)
    if(countpanel/=0) then
      call MPI_RECV(newpanel(1,mypanel),3*countpanel,MPI_INTEGER,nearproc(neighbor),0,comm3d,status,ierr)
    end if
  end do
end do
! curvatures
do neighbor=2,26,2
  do countobj=nobj4send(neighbor-1)+1,nobj4send(neighbor)
    myobj=locindex4sendobj(countobj)
    mypanel=npanel4proc(myobj-1)+1
    countpanel=npanel4proc(myobj)-npanel4proc(myobj-1)
    if(countpanel/=0) then
      call MPI_SEND(curvature(1,mypanel),3*countpanel,MPI_DOUBLE_PRECISION,nearproc(neighbor),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4recv(neighbor-2)+1,nobj4recv(neighbor-1)
    mypanel=npanel4recv(countobj-1)+1
    countpanel=npanel4recv(countobj)-npanel4recv(countobj-1)
    if(countpanel/=0) then
      call MPI_RECV(newcurvature(1,mypanel),3*countpanel,MPI_DOUBLE_PRECISION,nearproc(neighbor-1),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4send(neighbor-2)+1,nobj4send(neighbor-1)
    myobj=locindex4sendobj(countobj)
    mypanel=npanel4proc(myobj-1)+1
    countpanel=npanel4proc(myobj)-npanel4proc(myobj-1)
    if(countpanel/=0) then
      call MPI_SEND(newcurvature(1,mypanel),3*countpanel,MPI_DOUBLE_PRECISION,nearproc(neighbor-1),0,comm3d,status,ierr)
    end if
  end do
  do countobj=nobj4recv(neighbor-1)+1,nobj4recv(neighbor)
    mypanel=npanel4recv(countobj-1)+1
    countpanel=npanel4recv(countobj)-npanel4recv(countobj-1)
    if(countpanel/=0) then
      call MPI_RECV(newcurvature(1,mypanel),3*countpanel,MPI_DOUBLE_PRECISION,nearproc(neighbor),0,comm3d,status,ierr)
    end if
  end do
end do

deallocate(locindex4sendobj)

! save old object data
allocate(oldobj4proc(nobj4proc))
allocate(oldnvertex4proc(0:nobj4proc))
allocate(oldnpanel4proc(0:nobj4proc))
oldobj4proc=obj4proc
oldnvertex4proc=nvertex4proc
oldnpanel4proc=npanel4proc
deallocate(obj4proc)
deallocate(nvertex4proc)
deallocate(npanel4proc)

nlocvertex=oldnvertex4proc(nobj4proc)
nlocpanel=oldnpanel4proc(nobj4proc)

allocate(oldvertex(3,nlocvertex))
allocate(oldpanel(3,nlocpanel))
allocate(oldcurvature(3,nlocpanel))
oldvertex=vertex
oldpanel=panel
oldcurvature=curvature
deallocate(vertex)
deallocate(panel)
deallocate(curvature)

! count number of received repeated objects
nobj4repeat=0
do myobj=1,ntotalobj4recv
  do countobj=1,myobj-1
    if(obj4recv(1,countobj)/=0.and.obj4recv(1,myobj)==obj4recv(1,countobj)) then
      obj4recv(1,myobj)=0
      nobj4repeat=nobj4repeat+1
      exit
    end if
  end do
end do

! update local object info
oldnobj4proc=nobj4proc
nobj4proc=oldnobj4proc-nobj4delete+ntotalobj4recv-nobj4repeat

allocate(obj4proc(nobj4proc))
allocate(nvertex4proc(0:nobj4proc))
allocate(npanel4proc(0:nobj4proc))
obj4proc=0
nvertex4proc=0
npanel4proc=0

countobj=1
do myobj=1,oldnobj4proc
  if(proc4obj(0,myobj)==-1) then
    oldobj4proc(myobj)=0
  else
    obj4proc(countobj)=oldobj4proc(myobj)    
    nvertex4proc(countobj)=nvertex4proc(countobj-1)+(oldnvertex4proc(myobj)-oldnvertex4proc(myobj-1))
    npanel4proc(countobj)=npanel4proc(countobj-1)+(oldnpanel4proc(myobj)-oldnpanel4proc(myobj-1))
    countobj=countobj+1
  end if
end do
do myobj=1,ntotalobj4recv
  if(obj4recv(1,myobj)/=0) then
    obj4proc(countobj)=obj4recv(1,myobj)
    nvertex4proc(countobj)=nvertex4proc(countobj-1)+obj4recv(2,myobj)
    npanel4proc(countobj)=npanel4proc(countobj-1)+obj4recv(3,myobj)
    countobj=countobj+1
  end if
end do

deallocate(oldobj4proc)

nlocvertex=nvertex4proc(nobj4proc)
nlocpanel=npanel4proc(nobj4proc)

! update local vertices, panels and curvatures
allocate(vertex(3,nlocvertex))
allocate(panel(3,nlocpanel))
allocate(curvature(3,nlocpanel))

countobj=1
do myobj=1,oldnobj4proc
  if(proc4obj(0,myobj)/=-1) then
       vertex(:,nvertex4proc(countobj-1)+1:nvertex4proc(countobj))=&
    oldvertex(:,oldnvertex4proc(myobj-1)+1:oldnvertex4proc(myobj))
       panel(:,npanel4proc(countobj-1)+1:npanel4proc(countobj))=&
    oldpanel(:,oldnpanel4proc(myobj-1)+1:oldnpanel4proc(myobj))
       curvature(:,npanel4proc(countobj-1)+1:npanel4proc(countobj))=&
    oldcurvature(:,oldnpanel4proc(myobj-1)+1:oldnpanel4proc(myobj))
    countobj=countobj+1
  end if
end do

deallocate(proc4obj)
deallocate(oldnvertex4proc)
deallocate(oldnpanel4proc)
deallocate(oldvertex)
deallocate(oldpanel)
deallocate(oldcurvature)

do myobj=1,ntotalobj4recv
  if(obj4recv(1,myobj)/=0) then
       vertex(:,nvertex4proc(countobj-1)+1:nvertex4proc(countobj))=&
    newvertex(:,nvertex4recv(myobj-1)+1:nvertex4recv(myobj))
       panel(:,npanel4proc(countobj-1)+1:npanel4proc(countobj))=&
    newpanel(:,npanel4recv(myobj-1)+1:npanel4recv(myobj))
       curvature(:,npanel4proc(countobj-1)+1:npanel4proc(countobj))=&
    newcurvature(:,npanel4recv(myobj-1)+1:npanel4recv(myobj))
    countobj=countobj+1
  end if
end do

deallocate(obj4recv)
deallocate(nvertex4recv)
deallocate(npanel4recv)
deallocate(newvertex)
deallocate(newpanel)
deallocate(newcurvature)

! determine neighboring processors for each object
allocate(proc4obj(0:26,nobj4proc))
proc4obj=0

do myobj=1,nobj4proc
  do neighbor=0,26
    myproc=nearproc(neighbor)
    do mypanel=npanel4proc(myobj-1)+1,npanel4proc(myobj)
      A=vertex(:,nvertex4proc(myobj-1)+panel(1,mypanel))
      B=vertex(:,nvertex4proc(myobj-1)+panel(2,mypanel))
      C=vertex(:,nvertex4proc(myobj-1)+panel(3,mypanel))
      corner0=allcorner0(:,myproc)
      corner1=allcorner1(:,myproc)
      if(panelinPE(A,B,C,corner0,corner1)) then
        proc4obj(neighbor,myobj)=1
        exit ! regard the object on myproc if one of its panels on myproc
      endif
    end do
  end do
end do

end subroutine exchangeobject
