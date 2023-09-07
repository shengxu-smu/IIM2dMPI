subroutine locateobject
use MPI
use para
use Euler
use Lagrange
integer:: status(MPI_STATUS_SIZE)
integer:: i,countproc,countobj
integer:: myproc,myobj,myvertex,mypanel
integer:: nvertex,ndimension,npanel,npoint
integer:: nlocobj4myproc

integer, dimension(:), allocatable:: object
integer, dimension(:,:), allocatable:: processor4object
integer, dimension(:,:), allocatable:: object4processor
integer, dimension(:,:), allocatable:: panels

integer, dimension(:), allocatable:: nlocobj4proc
integer, dimension(:), allocatable:: objsend2proc

double precision, dimension(1:3):: A,B,C,corner0,corner1
double precision, dimension(:,:), allocatable:: vertices

character(5) filesuffix
namelist /vertexinfo/nvertex,ndimension
namelist /panelinfo/npanel,npoint

interface
  function panelinPE(A,B,C,corner0,corner1)
  double precision, dimension(1:3), intent(in):: A,B,C,corner0,corner1
  logical panelinPE
  end function panelinPE
end interface

allocate(object4processor(1:nobj,0:nprocs-1))
object4processor=0

if(myid==0) then
  allocate(object(nobj))
  allocate(processor4object(0:nprocs-1,1:nobj))
  object=0
  processor4object=-1

  open(100,file='INPUT/OBJECT/objlist')
  do myobj=1,nobj
    read(100,*)object(myobj)
    write(unit=filesuffix,fmt='(i0)') object(myobj)
    open(200,file='INPUT/OBJECT/vobj'//filesuffix)
    read(200,vertexinfo)
    allocate(vertices(3,nvertex))
    read(200,*)((vertices(i,myvertex),i=1,3),myvertex=1,nvertex)
    close(200)
    open(300,file='INPUT/OBJECT/pobj'//filesuffix)
    read(300,panelinfo)
    allocate(panels(3,npanel))
    read(300,*)((panels(i,mypanel),i=1,3),mypanel=1,npanel)
    close(300)

    countproc=0 ! count number of processors for object
    do myproc=0,nprocs-1
      do mypanel=1,npanel
        A=vertices(:,panels(1,mypanel))
        B=vertices(:,panels(2,mypanel))
        C=vertices(:,panels(3,mypanel))
        corner0=allcorner0(:,myproc)
        corner1=allcorner1(:,myproc)
        if(panelinPE(A,B,C,corner0,corner1)) then
          processor4object(countproc,myobj)=myproc
          countproc=countproc+1
          exit ! regard object on myproc if one of its panels on myproc
        endif
      end do
    end do
    deallocate(vertices)
    deallocate(panels)
  end do
  close(100)

  do myobj=1,nobj
    do countproc=0,nprocs-1
      if(processor4object(countproc,myobj)/=-1) then
        myproc=processor4object(countproc,myobj)
        do countobj=1,nobj
          if(object4processor(countobj,myproc)==0) then
            object4processor(countobj,myproc)=object(myobj)
            exit ! push object into stack
          end if
        end do
      end if
    end do
  end do
  deallocate(object)
  deallocate(processor4object)
end if

allocate(nlocobj4proc(0:nprocs-1))
if(myid==0) then
  do myproc=0,nprocs-1
    nlocobj4proc(myproc)=0
    do myobj=1,nobj  
     if(object4processor(myobj,myproc)/=0) then
        nlocobj4proc(myproc)=nlocobj4proc(myproc)+1
      end if
    end do
  end do
end if
call MPI_SCATTER(nlocobj4proc,1,MPI_INTEGER,nlocobj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

if(myid==0) then
  do myproc=1,nprocs-1
    nlocobj4myproc=nlocobj4proc(myproc)
    allocate(objsend2proc(nlocobj4myproc))
    objsend2proc=object4processor(1:nlocobj4myproc,myproc)
    call MPI_SEND(objsend2proc,nlocobj4myproc,MPI_INTEGER,myproc,0,MPI_COMM_WORLD,ierr)
    deallocate(objsend2proc)
  end do
  allocate(obj4proc(nlocobj))
  obj4proc=object4processor(1:nlocobj,0)
  deallocate(object4processor)
else
  allocate(obj4proc(nlocobj))
  call MPI_RECV(obj4proc,nlocobj,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
end if
deallocate(nlocobj4proc)

!---memory-inefficient option with fixed array size
!allocate(obj4proc(nobj))
!obj4proc=0
!call MPI_SCATTER(object4processor,nobj,MPI_INTEGER,obj4proc,nobj,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!deallocate(object4processor)

end subroutine locateobject
