SUBROUTINE locateobject
  USE MPI
  USE para
  USE Euler
  USE Lagrange
  INTEGER:: status(MPI_STATUS_SIZE)
  INTEGER:: countproc,countobj
  INTEGER:: myproc,myobj,myvertex,mypanel
  INTEGER:: nvertex,ndimension,npanel,npoint
  INTEGER:: nlocobj4myproc

  INTEGER, DIMENSION(:), ALLOCATABLE:: object
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: processor4object
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: object4processor
  INTEGER, DIMENSION(:,:), ALLOCATABLE:: panels

  INTEGER, DIMENSION(:), ALLOCATABLE:: nobjs4procs
  INTEGER, DIMENSION(:), ALLOCATABLE:: objsend2proc

  DOUBLE PRECISION, DIMENSION(1:2):: A,B,corner0,corner1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: vertices

  CHARACTER(5) filesuffix
  NAMELIST /objectinfo/nobj
  NAMELIST /vertexinfo/nvertex,ndimension
  NAMELIST /panelinfo/npanel,npoint

  INTERFACE
     FUNCTION panelinPE(A,B,corner0,corner1)
       DOUBLE PRECISION, DIMENSION(1:2), INTENT(in):: A,B,corner0,corner1
       LOGICAL panelinPE
     END FUNCTION panelinPE
  END INTERFACE

  OPEN(100,file='INPUT/OBJECT/objlist')
  READ(100,objectinfo)
  ALLOCATE(object(nobj))
  READ(100,*)(object(myobj),myobj=1,nobj)
  CLOSE(100)

  ALLOCATE(object_list(nobj))
  object_list = object

  ALLOCATE(object4processor(1:nobj,0:nprocs-1))
  object4processor=0
  IF(myid==0) THEN
     ALLOCATE(processor4object(0:nprocs-1,1:nobj))
     processor4object=-1
     DO myobj=1,nobj
        WRITE(unit=filesuffix,fmt='(i0)') object(myobj)
        OPEN(200,file='INPUT/OBJECT/vobj'//filesuffix)
        READ(200,vertexinfo)
        ALLOCATE(vertices(3,0:nvertex))
        READ(200,*)((vertices(i,myvertex),i=1,3),myvertex=0,nvertex)
        CLOSE(200)

        countproc=0 ! count number of processors for object
        DO myproc=0,nprocs-1
           DO myvertex=0,nvertex-1
              A=vertices(1:2,myvertex)
              B=vertices(1:2,myvertex+1)
              corner0=allcorner0(:,myproc)
              corner1=allcorner1(:,myproc)
              !print*, myproc,myvertex,A(1),B(1),corner0(1),corner1(1)
              IF(panelinPE(A,B,corner0,corner1)) THEN
                 processor4object(countproc,myobj)=myproc
                 countproc=countproc+1
                 EXIT ! regard object on myproc if one of its panels on myproc
              ENDIF
           END DO
        END DO
        DEALLOCATE(vertices)
     END DO

     DO myobj=1,nobj
        DO countproc=0,nprocs-1
           IF(processor4object(countproc,myobj)/=-1) THEN
              myproc=processor4object(countproc,myobj)
              DO countobj=1,nobj
                 IF(object4processor(countobj,myproc)==0) THEN
                    object4processor(countobj,myproc)=object(myobj)
                    EXIT ! push object into stack
                 END IF
              END DO
           END IF
        END DO
     END DO
     DEALLOCATE(object)
     DEALLOCATE(processor4object)
  END IF

  ALLOCATE(nobjs4procs(0:nprocs-1))
  IF(myid==0) THEN
     DO myproc=0,nprocs-1
        nobjs4procs(myproc)=0
        DO myobj=1,nobj
           IF(object4processor(myobj,myproc)/=0) THEN
              nobjs4procs(myproc)=nobjs4procs(myproc)+1
           END IF
        END DO
     END DO
  END IF
  CALL MPI_SCATTER(nobjs4procs,1,MPI_INTEGER,nobj4proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  IF(myid==0) THEN
     DO myproc=1,nprocs-1
        nlocobj4myproc=nobjs4procs(myproc)
        ALLOCATE(objsend2proc(nlocobj4myproc))
        objsend2proc=object4processor(1:nlocobj4myproc,myproc)
        CALL MPI_SEND(objsend2proc,nlocobj4myproc,MPI_INTEGER,myproc,0,MPI_COMM_WORLD,ierr)
        DEALLOCATE(objsend2proc)
     END DO
     ALLOCATE(obj4proc(nobj4proc))
     obj4proc=object4processor(1:nobj4proc,0)
     DEALLOCATE(object4processor)
  ELSE
     ALLOCATE(obj4proc(nobj4proc))
     CALL MPI_RECV(obj4proc,nobj4proc,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
  END IF
  DEALLOCATE(nobjs4procs)


END SUBROUTINE locateobject
