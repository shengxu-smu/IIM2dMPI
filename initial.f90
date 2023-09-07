SUBROUTINE initial

  USE para
  USE Euler
  USE Lagrange

  object_move =.FALSE.
  IF(isingular.EQV..TRUE.) THEN
     ! allocate moving object properties

     ! allocate moving properties
     ALLOCATE(xsc(1:npbj4proc))
     ALLOCATE(ysc(1:nobj4proc))
     ALLOCATE(theta(1:nobj4proc))
     ALLOCATE(xsct(1:npbj4proc))
     ALLOCATE(ysct(1:nobj4proc))
     ALLOCATE(thetat(1:nobj4proc))
     ALLOCATE(xsctt(1:npbj4proc))
     ALLOCATE(ysctt(1:nobj4proc))
     ALLOCATE(thetatt(1:nobj4proc))

     ALLOCATE(us(1:nvertex4proc(nobj4proc)))
     ALLOCATE(vs(1:nvertex4proc(nobj4proc)))
     ALLOCATE(taox(1:nvertex4proc(nobj4proc)))
     ALLOCATE(taoy(1:nvertex4proc(nobj4proc)))
  ENDIF

  nstart=0
  t=ZERO

  !dt0 =0.5d0
  !cflv=0.5d0
  !cflc=0.5d0

  ! setup initial movement of object
  ! velocity

  ! Setup initial u,v,p,d,dn
  ! Setup BCs
  ! Setup uBC
  ! Setup vBC
  ! Setup pBC
  ! exchange ghost point

  !iop=0
  !iou=0
  !iov=0

  p=ZERO
  d=ZERO
  dn=ZERO
  u=ZERO
  v=ZERO
  !u=dble(coords(1)+0.1*coords(2))
  !v=dble(coords(2)+0.2*coords(2))

  ! initial flow field
  DO j=jubeg,juend
     DO i=iubeg,iuend
        !    u(i,j)=2.0d0*(i+1)+0.1d0*(j+1)
        u(i,j)=1.0d0
     ENDDO
  ENDDO
  DO j=jvbeg,jvend
     DO i=ivbeg,ivend
        !    v(i,j)=-i-0.2d0*j
        v(i,j)=0.0d0
     ENDDO
  ENDDO

  CALL exchange4ghost

  ! initial all right hand side
  rhsp=ZERO
  rhsu1=ZERO
  rhsv1=ZERO
  rhsu2=ZERO
  rhsv2=ZERO
  rhsu3=ZERO
  rhsv3=ZERO
  rhsu4=ZERO
  rhsv4=ZERO

  ! initial velocity with jump conditions

  CALL velocity_initial

END SUBROUTINE initial
