SUBROUTINE jc_firstsecond

  USE mpi
  USE para
  USE stretchxy
  USE Euler
  USE Lagrange

  INTEGER:: myvertex,myobj

  DOUBLE PRECISION:: gacobi2,gacobi,tx,ty,xn,yn
  DOUBLE PRECISION:: duxjc,duyjc,dduxjc,dduyjc,dduxyjc,ddpxjc0,ddpxjc1
  DOUBLE PRECISION:: dvxjc,dvyjc,ddvxjc,ddvyjc,ddvxyjc,ddpyjc0,ddpyjc1
  DOUBLE PRECISION:: dpxjc,dpyjc
  DOUBLE PRECISION:: dudnjc,dudnjcn,dudnjcp,dvdnjc,dvdnjcn,dvdnjcp
  DOUBLE PRECISION:: ddudnp,ddvdnp,qx,qy,xs0,ys0,xs1,ys1

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE:: deltaujc1,deltavjc1,deltaujc2,deltavjc2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE:: ujc2,vjc2,pjc2

  ! dudxdt: derivative of dudx over tao
  DOUBLE PRECISION:: dudxdt,dvdxdt,dpdxdt,dudydt,dvdydt,dpdydt

  ! variables used in interpolation of dudnjc ddudnjc
  DOUBLE PRECISION:: ds,dist,foo,r3
  DOUBLE PRECISION,DIMENSION(0:3):: uu,vv,xx,yy,distance
  INTEGER:: ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
  INTEGER, PARAMETER:: many=3,nany=3
  DOUBLE PRECISION,DIMENSION(1:many):: xa,ya,xb,yb
  DOUBLE PRECISION, DIMENSION(1:2):: A,B,corner0,corner1

  ! pressure solver test
  DOUBLE PRECISION:: dpdnp,dpdnn,dpdnjc,dpdxtrue,dpdytrue,d2pdx2true,d2pdp2true
  DOUBLE PRECISION:: jcp1,jcp2,dpdt
  DOUBLE PRECISION:: err1,err2, maxerr1,maxerr2,maxerr

  maxerr1 = -100.0d0
  maxerr2 = -100.0d0

  ujc = zero
  vjc = zero
  pjc = zero

  ALLOCATE(deltaujc1(1:nvertex4proc(nobj4proc)))
  ALLOCATE(deltavjc1(1:nvertex4proc(nobj4proc)))
  ALLOCATE(deltaujc2(1:nvertex4proc(nobj4proc)))
  ALLOCATE(deltavjc2(1:nvertex4proc(nobj4proc)))

  ALLOCATE(ujc2(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(vjc2(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(pjc2(1:8, 1:nvertex4proc(nobj4proc)))

  us = zero
  vs = zero
  xx = zero
  yy = zero
  uu = zero
  vv = zero
  distance = zero
  thetat = zero
  thetatt = zero

  ds=1.01d0*SQRT(dxi*dxi+deta*deta)
  ! counter clockwise direction
  DO myobj=1,nobj4proc
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
        A = vertex(:,myvertex)
        IF(myvertex.NE.nvertex4proc(myobj)) THEN
           B = vertex(:,myvertex+1)
        ELSE
           B = vertex(:,nvertex4proc(myobj-1)+2)
        ENDIF
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           ! tao direction
           tx=taox(myvertex)
           ty=taoy(myvertex)

           gacobi=dsqrt(tx*tx+ty*ty)
           gacobi2=tx*tx+ty*ty

           ! normal direction
           xn= ty
           yn=-tx
           xn=xn/gacobi
           yn=yn/gacobi

           xx(0)=vertex(1,myvertex)
           yy(0)=vertex(2,myvertex)

           jcp1 = dsin(A(1))*dsin(A(2))-dexp(-A(1))*dexp(-A(2))
           jcp2 = dsin(B(1))*dsin(B(2))-dexp(-B(1))*dexp(-B(2))
           dist = SQRT((B(1)-A(1))*(B(1)-A(1))+(B(2)-A(2))*(B(2)-A(2)))
           dpdt = (jcp2-jcp1)/dist

           dpdt = dcos(xx(0))*dsin(yy(0))*tx+dsin(xx(0))*dcos(yy(0))*ty + &
                  dexp(-xx(0))*dexp(-yy(0))*(tx+ty)

           dpdnp = dcos(xx(0))*dsin(yy(0))*xn+dsin(xx(0))*dcos(yy(0))*yn
           dpdnn = -dexp(-xx(0))*dexp(-yy(0))*(xn+yn)
           dpdnjc = dpdnp-dpdnn
           dpxjc = (yn*dpdt-ty*dpdnjc)/(tx*yn-ty*xn)
           dpyjc = (tx*dpdnjc-xn*dpdt)/(tx*yn-ty*xn)

           pjc(1,myvertex)=dpxjc
           pjc(2,myvertex)=dpyjc

           !  print*,myid,myvertex,dpxjc,dpyjc
        ENDIF
     ENDDO ! end myvertex

     DO myvertex=nvertex4proc(myobj),nvertex4proc(myobj-1)+1,-1
        A = vertex(:,myvertex)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN
           ! tao direction
           IF ( myvertex ==  nvertex4proc(myobj-1)+1) THEN
              B = vertex(:,nvertex4proc(myobj)-1)
              tx = -taox(nvertex4proc(myobj)-1)
              ty = -taoy(nvertex4proc(myobj)-1)
           ELSE
              B = vertex(:,myvertex-1)
              tx = -taox(myvertex-1)
              ty = -taoy(myvertex-1)
           END IF

           gacobi=SQRT(tx*tx+ty*ty)
           gacobi2=tx*tx+ty*ty

           ! normal direction
           xn=-ty
           yn= tx
           xn=xn/gacobi
           yn=yn/gacobi

           xx(0)=vertex(1,myvertex)
           yy(0)=vertex(2,myvertex)

           jcp1 = dsin(A(1))*dsin(A(2))-dexp(-A(1))*dexp(-A(2))
           jcp2 = dsin(B(1))*dsin(B(2))-dexp(-B(1))*dexp(-B(2))
           dist = SQRT((B(1)-A(1))*(B(1)-A(1))+(B(2)-A(2))*(B(2)-A(2)))
           dpdt = (jcp2-jcp1)/dist

           dpdt = dcos(xx(0))*dsin(yy(0))*tx+dsin(xx(0))*dcos(yy(0))*ty + &
                  dexp(-xx(0))*dexp(-yy(0))*(tx+ty)

           dpdnp = dcos(xx(0))*dsin(yy(0))*xn+dsin(xx(0))*dcos(yy(0))*yn
           dpdnn = -dexp(-xx(0))*dexp(-yy(0))*(xn+yn)
           dpdnjc = dpdnp-dpdnn
           dpxjc = (yn*dpdt-ty*dpdnjc)/(tx*yn-ty*xn)
           dpyjc = (tx*dpdnjc-xn*dpdt)/(tx*yn-ty*xn)

           pjc2(1,myvertex)=dpxjc
           pjc2(2,myvertex)=dpyjc
        ENDIF
     ENDDO

     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)

        A = vertex(:,myvertex)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           xx(0)=vertex(1,myvertex)
           yy(0)=vertex(2,myvertex)
           pjc(1,myvertex)=(pjc(1,myvertex)+pjc2(1,myvertex))/2.0d0
           pjc(2,myvertex)=(pjc(2,myvertex)+pjc2(2,myvertex))/2.0d0

           dpdxtrue = dcos(xx(0))*dsin(yy(0))+dexp(-xx(0))*dexp(-yy(0))
           dpdytrue = dsin(xx(0))*dcos(yy(0))+dexp(-xx(0))*dexp(-yy(0))

           print*,myvertex,pjc(2,myvertex),dsin(xx(0))*dcos(yy(0))+dexp(-xx(0))*dexp(-yy(0)),&
                    ABS(pjc(2,myvertex)-dpdytrue)
           maxerr1 = MAX(maxerr1,ABS(pjc(1,myvertex)-dpdxtrue))
           maxerr2 = MAX(maxerr2,ABS(pjc(2,myvertex)-dpdytrue))

        ENDIF
     ENDDO
  ENDDO
  PRINT*,'error of dpdx,dpdy',myid,maxerr1,maxerr2


  DO myobj=1, nobj4proc

     ! second order jump conditions
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
        A=vertex(:,myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           ! tao direction
           tx=taox(myvertex)
           ty=taoy(myvertex)
           gacobi2=tx*tx+ty*ty

           ! dudxdt, dvdxdt, dudydt, dxdydt
           ! dudxdt = (dudxB-dudxA)/|AB|
           ! position of neighbor vertex
           xs0=vertex(1,myvertex)
           xs1=vertex(1,myvertex+1)
           ys0=vertex(2,myvertex)
           ys1=vertex(2,myvertex+1)
           dist=SQRT((ys1-ys0)*(ys1-ys0)+(xs1-xs0)*(xs1-xs0))

           dpdxdt = (pjc(1,myvertex+1) - pjc(1,myvertex))/dist
           dpdydt = (pjc(2,myvertex+1) - pjc(2,myvertex))/dist

           !  print*,myid,myvertex,dpdxdt,pjc(1,myvertex+1),pjc(1,myvertex)
           ! ddpdx ddpdy jump conditions
           r3 = -2.0d0*dsin(A(1))*dsin(A(2))-2.0d0*dexp(-A(1))*dexp(-A(2))
           ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2
           ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           !  print*,myid,myvertex,ddpxjc0,ddpyjc0,dpdxdt,dpdydty,r3
          !  print*,myvertex,pjc(2,myvertex+1),pjc(2,myvertex),dpdydt,dist

           pjc(3,myvertex)=ddpxjc1
           pjc(7,myvertex)=ddpxjc0
           pjc(4,myvertex)=ddpyjc1
           pjc(8,myvertex)=ddpyjc0

           dpdxtrue = -dsin(A(1))*dsin(A(2))-dexp(-A(1))*dexp(-A(2))
           dpdytrue = -dsin(A(1))*dsin(A(2))-dexp(-A(1))*dexp(-A(2))

          !  print*,myvertex,pjc(3,myvertex),dpdxtrue
        ENDIF
     ENDDO !end my vertex
     pjc(3,nvertex4proc(myobj)) = pjc(3,nvertex4proc(myobj-1)+1)
     pjc(7,nvertex4proc(myobj)) = pjc(7,nvertex4proc(myobj-1)+1)
     pjc(4,nvertex4proc(myobj)) = pjc(4,nvertex4proc(myobj-1)+1)
     pjc(8,nvertex4proc(myobj)) = pjc(8,nvertex4proc(myobj-1)+1)

     ! second order jump conditions
     DO myvertex=nvertex4proc(myobj),nvertex4proc(myobj-1)+2,-1
        A=vertex(:,myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           ! tao direction
           tx = -taox(myvertex-1)
           ty = -taoy(myvertex-1)
           gacobi2=tx*tx+ty*ty

           ! dudxdt, dvdxdt, dudydt, dxdydt
           ! dudxdt = (dudxB-dudxA)/|AB|
           ! position of neighbor vertex
           xs0=vertex(1,myvertex)
           xs1=vertex(1,myvertex-1)
           ys0=vertex(2,myvertex)
           ys1=vertex(2,myvertex-1)
           dist=SQRT((ys1-ys0)*(ys1-ys0)+(xs1-xs0)*(xs1-xs0))

           dpdxdt = (pjc2(1,myvertex-1) - pjc2(1,myvertex))/dist
           dpdydt = (pjc2(2,myvertex-1) - pjc2(2,myvertex))/dist

           ! ddpdx ddpdy jump conditions
           r3 = -2.0d0*dsin(A(1))*dsin(A(2))-2.0d0*dexp(-A(1))*dexp(-A(2))
           ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2
           ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           pjc2(3,myvertex)=ddpxjc1
           pjc2(7,myvertex)=ddpxjc0
           pjc2(4,myvertex)=ddpyjc1
           pjc2(8,myvertex)=ddpyjc0

           !  print*,myid,myvertex,pjc2(3,myvertex),pjc2(4,myvertex)

        ENDIF
     ENDDO ! end myvertex
     pjc2(3,nvertex4proc(myobj-1)+1) = pjc2(3,nvertex4proc(myobj))
     pjc2(7,nvertex4proc(myobj-1)+1) = pjc2(7,nvertex4proc(myobj))
     pjc2(4,nvertex4proc(myobj-1)+1) = pjc2(4,nvertex4proc(myobj))
     pjc2(8,nvertex4proc(myobj-1)+1) = pjc2(8,nvertex4proc(myobj))

     maxerr1=-100
     maxerr2=-100
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)

        A=vertex(:,myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)
        corner0=(/xf(iubeg-1),yf(jvbeg-1)/)
        corner1=(/xf(iuend+1),yf(jvend+1)/)

        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN

           xx(0)=vertex(1,myvertex)
           yy(0)=vertex(2,myvertex)

           !  print*,myid,myvertex,pjc(7,myvertex),pjc2(7,myvertex)


           pjc(3,myvertex)=(pjc(3,myvertex)+pjc2(3,myvertex))/2.0d0
           pjc(4,myvertex)=(pjc(4,myvertex)+pjc2(4,myvertex))/2.0d0
           pjc(7,myvertex)=(pjc(7,myvertex)+pjc2(7,myvertex))/2.0d0
           pjc(8,myvertex)=(pjc(8,myvertex)+pjc2(8,myvertex))/2.0d0

           dpdxtrue = -dsin(A(1))*dsin(A(2))-dexp(-A(1))*dexp(-A(2))
           dpdytrue = -dsin(A(1))*dsin(A(2))-dexp(-A(1))*dexp(-A(2))

           maxerr1 = MAX(maxerr1,ABS(pjc(7,myvertex)-dpdxtrue))
           maxerr2 = MAX(maxerr2,ABS(pjc(8,myvertex)-dpdytrue))
          !  print*,myvertex,pjc(7,myvertex),dpdxtrue,pjc(8,myvertex),dpdytrue

        ENDIF
     ENDDO
  ENDDO
  PRINT*,'error of ddpdx,ddpdy',myid,maxerr1,maxerr2

  DEALLOCATE(deltaujc1)
  DEALLOCATE(deltavjc1)
  DEALLOCATE(deltaujc2)
  DEALLOCATE(deltavjc2)

  DEALLOCATE(ujc2)
  DEALLOCATE(vjc2)
  DEALLOCATE(pjc2)

END SUBROUTINE jc_firstsecond
