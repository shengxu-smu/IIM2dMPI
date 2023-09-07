SUBROUTINE initial

  USE para
  USE Euler
  USE Lagrange

  ! allocate moving object properties
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

  !icjc=0
  !ifjc=0
  !icjf=0

  !ickc=0
  !ifkc=0
  !ickf=0

  !jckc=0
  !jfkc=0
  !jckf=0

  p=ZERO
  d=ZERO
  dn=ZERO
  u=ZERO
  v=ZERO
  !u=dble(coords(1)+0.1*coords(2))
  !v=dble(coords(2)+0.2*coords(2))
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
  !call interpolateuv
  !      write(*,*),'ucc',ucc(1,1),ucc(1,101),ucc(101,1),ucc(101,101)
  !      write(*,*),'vcc',vcc(50,50),'vec',vec(50,50),'vee',vee(50,50)
  !call pbc
  !print*,u(iubeg,jubeg),u(iuend,juend),v(ivbeg,jvbeg),v(ivend,jvend)

  rhsp=ZERO
  rhsu1=ZERO
  rhsv1=ZERO
  rhsu2=ZERO
  rhsv2=ZERO
  rhsu3=ZERO
  rhsv3=ZERO
  rhsu4=ZERO
  rhsv4=ZERO


  CALL velocity_initial


END SUBROUTINE initial
