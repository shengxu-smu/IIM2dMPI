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

  ! pressure solver test
  DOUBLE PRECISION:: dpdnp,dpdnn,dpdnjc,dpdxtrue,dpdytrue,d2pdx2true,d2pdp2true
  DOUBLE PRECISION:: err1,err2, maxerr1,maxerr2,maxerr

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
  DOUBLE PRECISION, DIMENSION(1:2):: A,corner0,corner1

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

           ! velocity and position of vertices
           uu(0)=us(myvertex)
           vv(0)=vs(myvertex)
           xx(0)=vertex(1,myvertex)
           yy(0)=vertex(2,myvertex)
           distance(0)=zero

           dpdnp = dcos(xx(0))*dsin(yy(0))*xn+dsin(xx(0))*dcos(yy(0))*yn
           dpdnn = -dexp(-xx(0))*dexp(-yy(0))*(xn+yn)
           dpdnjc = dpdnp-dpdnn
           dpxjc = -ty/(tx*yn-ty*xn)*dpdnjc
           dpyjc =  tx/(tx*yn-ty*xn)*dpdnjc

           dpdxtrue = dcos(xx(0))*dsin(yy(0))+dexp(-xx(0))*dexp(-yy(0))
           dpdytrue = dsin(xx(0))*dcos(yy(0))+dexp(-xx(0))*dexp(-yy(0))

           PRINT*,myvertex,dpxjc,dpdxtrue,dpyjc,dpdytrue

           maxerr1 = MAX(maxerr1,dpxjc-dpdxtrue)
           maxerr2 = MAX(maxerr2,dpyjc-dpdytrue)

           ujc(1,myvertex)=duxjc
           vjc(1,myvertex)=dvxjc
           pjc(1,myvertex)=dpxjc
           ujc(2,myvertex)=duyjc
           vjc(2,myvertex)=dvyjc
           pjc(2,myvertex)=dpyjc
           !  PRINT*, myid, myvertex,ujc(1,myvertex),vjc(1,myvertex),pjc(1,myvertex)
        ENDIF !
     ENDDO ! end myvertex
     PRINT*,maxerr2,maxerr1

     ! second order jump conditions
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
        A=vertex(:,myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)

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

           dudxdt = (ujc(1,myvertex+1) - ujc(1,myvertex))/dist
           dvdxdt = (vjc(1,myvertex+1) - vjc(1,myvertex))/dist
           dpdxdt = (pjc(1,myvertex+1) - pjc(1,myvertex))/dist

           dudydt = (ujc(2,myvertex+1) - ujc(2,myvertex))/dist
           dvdydt = (vjc(2,myvertex+1) - vjc(2,myvertex))/dist
           dpdydt = (pjc(2,myvertex+1) - pjc(2,myvertex))/dist


           ! dduxjc dduyjc dduxyjc
           dduxjc=(tx*dudxdt-ty*dudydt+ty*ty*deltaujc1(myvertex))/gacobi2
           ddvxjc=(tx*dvdxdt-ty*dvdydt+ty*ty*deltavjc1(myvertex))/gacobi2

           dduyjc=(ty*dudydt-tx*dudxdt+tx*tx*deltaujc1(myvertex))/gacobi2
           ddvyjc=(ty*dvdydt-tx*dvdxdt+tx*tx*deltavjc1(myvertex))/gacobi2

           dduxyjc=(ty*dudxdt+tx*dudydt-tx*ty*deltaujc1(myvertex))/gacobi2
           ddvxyjc=(ty*dvdxdt+tx*dvdydt-tx*ty*deltavjc1(myvertex))/gacobi2

           ! ddpdx ddpdy jump conditions
           r3 = -2.0d0*dsin(xs0)*dsin(ys0)-2.0d0*dexp(-xs0)*dexp(-ys0)
           !  r3 = 0.0d0
           ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           !  r3 = 1.0d0
           ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           ujc(3,myvertex)=dduxjc
           vjc(3,myvertex)=ddvxjc
           pjc(3,myvertex)=ddpxjc1
           ujc(5,myvertex)=dduxyjc
           vjc(5,myvertex)=ddvxyjc
           pjc(7,myvertex)=ddpxjc0

           ujc(4,myvertex)=dduyjc
           vjc(4,myvertex)=ddvyjc
           pjc(4,myvertex)=ddpyjc1
           ujc(6,myvertex)=dduxyjc
           vjc(6,myvertex)=ddvxyjc
           pjc(8,myvertex)=ddpyjc0
        ENDIF
     ENDDO !end my vertex
     ujc(3,nvertex4proc(myobj)) = ujc(3,nvertex4proc(myobj-1)+1)
     vjc(3,nvertex4proc(myobj)) = vjc(3,nvertex4proc(myobj-1)+1)
     pjc(3,nvertex4proc(myobj)) = pjc(3,nvertex4proc(myobj-1)+1)
     ujc(5,nvertex4proc(myobj)) = ujc(5,nvertex4proc(myobj-1)+1)
     vjc(5,nvertex4proc(myobj)) = vjc(5,nvertex4proc(myobj-1)+1)
     pjc(7,nvertex4proc(myobj)) = pjc(7,nvertex4proc(myobj-1)+1)

     ujc(4,nvertex4proc(myobj)) = ujc(4,nvertex4proc(myobj-1)+1)
     vjc(4,nvertex4proc(myobj)) = vjc(4,nvertex4proc(myobj-1)+1)
     pjc(4,nvertex4proc(myobj)) = pjc(4,nvertex4proc(myobj-1)+1)
     ujc(6,nvertex4proc(myobj)) = ujc(6,nvertex4proc(myobj-1)+1)
     vjc(6,nvertex4proc(myobj)) = vjc(6,nvertex4proc(myobj-1)+1)
     pjc(8,nvertex4proc(myobj)) = pjc(8,nvertex4proc(myobj-1)+1)

  ENDDO ! end myobj

  DEALLOCATE(deltaujc1)
  DEALLOCATE(deltavjc1)
  DEALLOCATE(deltaujc2)
  DEALLOCATE(deltavjc2)

  DEALLOCATE(ujc2)
  DEALLOCATE(vjc2)
  DEALLOCATE(pjc2)



END SUBROUTINE jc_firstsecond
