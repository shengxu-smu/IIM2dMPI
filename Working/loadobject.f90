subroutine loadobject
use para
use Euler
use Lagrange
integer:: i,countvertex,countpanel,countproc
integer:: myobj,myvertex,mypanel,myproc
integer:: object,nvertex,ndimension,npanel,npoint

double precision, dimension(1:3):: A,B,C,corner0,corner1

character(5) filesuffix
namelist /vertexinfo/nvertex,ndimension
namelist /panelinfo/npanel,npoint

interface
  function panelinPE(A,B,C,corner0,corner1)
  double precision, dimension(1:3), intent(in):: A,B,C,corner0,corner1
  logical panelinPE
  end function panelinPE
end interface

!nlocobj=0
!do myobj=1,nobj  
!  if(obj4proc(myobj)/=0) then
!    nlocobj=nlocobj+1
!  end if
!end do

allocate(nvertex4obj(nlocobj))
allocate(npanel4obj(nlocobj))

nlocvertex=0
nlocpanel=0
do myobj=1,nlocobj
  object=obj4proc(myobj)
  write(filesuffix,'(i0)') object
  open(100,file='INPUT/OBJECT/vobj'//filesuffix)
  read(100,vertexinfo)
  nvertex4obj(myobj)=nvertex
  nlocvertex=nlocvertex+nvertex
  close(100)
  open(200,file='INPUT/OBJECT/pobj'//filesuffix)
  read(200,panelinfo)
  npanel4obj(myobj)=npanel
  nlocpanel=nlocpanel+npanel
  close(200)
end do
allocate(vertex(3,nlocvertex))
allocate(panel(3,nlocpanel))
allocate(curvature(3,nlocpanel))

countvertex=0
countpanel=0
do myobj=1,nlocobj
  object=obj4proc(myobj)
  write(filesuffix,'(i0)') object
  open(100,file='INPUT/OBJECT/vobj'//filesuffix)
  read(100,vertexinfo)
  read(100,*)((vertex(i,countvertex+myvertex),i=1,3),myvertex=1,nvertex)
  countvertex=countvertex+nvertex
  close(100)
  open(200,file='INPUT/OBJECT/pobj'//filesuffix)
! open(300,file='INPUT/OBJECT/cobj'//filesuffix)
  read(200,panelinfo)
! read(300,panelinfo)
  read(200,*)((panel(i,countpanel+mypanel),i=1,3),mypanel=1,npanel)
! read(300,*)((allcurvature(i,countpanel+mypanel),i=1,3),mypanel=1,npanel)
  countpanel=countpanel+npanel
  close(200)
! close(300)
end do

allocate(nearproc4obj(26,nlocobj))

countvertex=0
countpanel=0
nlocproc=0
do myobj=1,nlocobj
  nvertex=nvertex4obj(myobj)
  npanel=npanel4obj(myobj)
  do myproc=0,nprocs-1
    do mypanel=countpanel+1,countpanel+npanel
      A=vertex(:,countvertex+panel(1,mypanel))
      B=vertex(:,countvertex+panel(2,mypanel))
      C=vertex(:,countvertex+panel(3,mypanel))
      corner0=allcorner0(:,myproc)
      corner1=allcorner1(:,myproc)
      if(panelinPE(A,B,C,corner0,corner1)) then
        nlocproc=nlocproc+1
        exit ! regard the object on myproc if one of its panels on myproc
      endif
    end do
  end do
  countvertex=countvertex+nvertex
  countpanel=countpanel+npanel
end do
allocate(proc4obj(nlocproc))

countvertex=0
countpanel=0
countproc=0
do myobj=1,nlocobj
  nvertex=nvertex4obj(myobj)
  npanel=npanel4obj(myobj)
  do myproc=0,nprocs-1
    do mypanel=countpanel+1,countpanel+npanel
      A=vertex(:,countvertex+panel(1,mypanel))
      B=vertex(:,countvertex+panel(2,mypanel))
      C=vertex(:,countvertex+panel(3,mypanel))
      corner0=allcorner0(:,myproc)
      corner1=allcorner1(:,myproc)
      if(panelinPE(A,B,C,corner0,corner1)) then
        countproc=countproc+1
        proc4obj(countproc)=myproc
        exit
      endif
    end do
  end do
  countvertex=countvertex+nvertex
  countpanel=countpanel+npanel
end do

print *,'myid = ',myid,' nlocobj = ',nlocobj,' nlocvertex = ',nlocvertex, 'nlocpanel = ',nlocpanel, &
        ' nlocproc = ',nlocproc,' proc4obj = ',proc4obj

end subroutine loadobject
