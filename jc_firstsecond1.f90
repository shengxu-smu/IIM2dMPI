SUBROUTINE jc_firstsecond

  USE mpi
  USE para
  USE stretchxy
  USE Euler
  USE Lagrange

  INTEGER:: myvertex,myobj
  INTEGER:: mylocalvertex,mylocalvertex1
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
  DOUBLE PRECISION, DIMENSION(1:2):: A,corner0,corner1


  ALLOCATE(deltaujc1(1:nvertexlocal(nobj4proc)))
  ALLOCATE(deltavjc1(1:nvertexlocal(nobj4proc)))
  ALLOCATE(deltaujc2(1:nvertexlocal(nobj4proc)))
  ALLOCATE(deltavjc2(1:nvertexlocal(nobj4proc)))

  ALLOCATE(ujc2(1:6, 1:nvertexlocal(nobj4proc)))
  ALLOCATE(vjc2(1:6, 1:nvertexlocal(nobj4proc)))
  ALLOCATE(pjc2(1:8, 1:nvertexlocal(nobj4proc)))

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
     DO mylocalvertex=nvertexlocal(myobj-1)+1,nvertexlocal(myobj)
        !get index of vertices
        myvertex=vertexindexlocal(mylocalvertex)

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

        ! dudn jump conditions on negative side
        ! dudnjcn= thetat X n
        dudnjcn = -thetat(myobj)*yn
        dvdnjcn =  thetat(myobj)*xn

        ! dudn jump conditions on positive side
        ! ============ put this part into a subroutine ============
        ! a. find 3 points on normal direction s1, s2, s3

        DO n=1,nany

           ! coordinates of interpolation point
           xx(n)= xi2x( x2xi( vertex(1,myvertex)) +DBLE(n)*ds*xn)
           yy(n)= eta2y(y2eta(vertex(2,myvertex)) +DBLE(n)*ds*yn)
           distance(n)=SQRT((xx(n)-xx(0))*(xx(n)-xx(0))+(yy(n)-yy(0))*(yy(n)-yy(0)))

           ! indices of interpolation point
           i=INT((x2xi(xx(n))-x2xi(x0))/dxi*2)
           j=INT((y2eta(yy(n))-y2eta(y0))/deta*2)
           IF(MOD(i,2)==0) THEN
              ie=i/2
              ic=ie
           ELSE
              ic=(i+1)/2
              ie=ic-1
           ENDIF
           IF(MOD(j,2)==0) THEN
              je=j/2
              jc=je
           ELSE
              jc=(j+1)/2
              je=jc-1
           ENDIF
           iu=ie
           ju=jc
           iv=ic
           jv=je

           id=INT(SIGN(1.d0,xn))
           jd=INT(SIGN(1.d0,yn))
           IF(id<0.0d0) THEN
              iu=iu+1
              iv=iv+1
           ENDIF
           IF(jd<0.0d0) THEN
              ju=ju+1
              jv=jv+1
           ENDIF

           ! interpolation of uu(n)
           DO j=1,many
              DO i=1,many
                 xa(i)=xf(iu+id*(i-1))
                 ya(i)=u(iu+id*(i-1),ju+jd*(j-1))
              ENDDO
              xb(j)=yc(ju+jd*(j-1))
              CALL interpolate(xa,ya,many,xx(n),yb(j),foo)
           ENDDO
           CALL interpolate(xb,yb,many,yy(n),uu(n),foo)

           ! interpolation of vv(n)
           DO i=1,many
              DO j=1,many
                 xa(j)=yf(jv+jd*(j-1))
                 ya(j)=v(iv+id*(i-1),jv+jd*(j-1))
              ENDDO
              xb(i)=xc(iv+id*(i-1))
              CALL interpolate(xa,ya,many,yy(n),yb(i),foo)
           ENDDO
           CALL interpolate(xb,yb,many,xx(n),vv(n),foo)
        ENDDO

        ! b. dudn jump conditions on positive side
        CALL FiniteDiff(distance,uu,dudnjcp,ddudnp)
        CALL FiniteDiff(distance,vv,dvdnjcp,ddvdnp)

        ! c. dudn jump conditions
        dudnjc = dudnjcp - dudnjcn
        dvdnjc = dvdnjcp - dvdnjcn

! if(myid == 2 .and. mylocalvertex<=4) then
!   print*, mylocalvertex, vertexindexlocal(mylocalvertex),uu(0),uu(1),uu(2),uu(3)
! endif

        ! laplacian u jump conditions
        ! ddudn jump conditions
        deltaujc1(mylocalvertex) = ddudnp + vertex(3,myvertex)*dudnjc
        deltavjc1(mylocalvertex) = ddvdnp + vertex(3,myvertex)*dvdnjc

        ! dudx dvdx dudy dvdy
        duxjc = -ty/(tx*yn-ty*xn)*dudnjc
        dvxjc = -ty/(tx*yn-ty*xn)*dvdnjc
        duyjc =  tx/(tx*yn-ty*xn)*dudnjc
        dvyjc =  tx/(tx*yn-ty*xn)*dvdnjc

        ! dpdx dpdy jump conditions
        ! body force qx, qy
        qx =  thetatt(myobj)*(vertex(2,myvertex)-ysc(myobj))
        qy = -thetatt(myobj)*(vertex(1,myvertex)-xsc(myobj))

        ! dpdx, dpdy jump conditions
        dpxjc = 1.0d0/Re*deltaujc1(mylocalvertex) + qx
        dpyjc = 1.0d0/re*deltavjc1(mylocalvertex) + qy

        ujc(1,mylocalvertex)=duxjc
        vjc(1,mylocalvertex)=dvxjc
        pjc(1,mylocalvertex)=dpxjc
        ujc(2,mylocalvertex)=duyjc
        vjc(2,mylocalvertex)=dvyjc
        pjc(2,mylocalvertex)=dpyjc

     ENDDO

     ! second order jump conditions
     DO mylocalvertex=nvertexlocal(myobj-1)+1,nvertexlocal(myobj)
        myvertex=vertexindexlocal(mylocalvertex)
        A=vertex(:,myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)

        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN
           IF(myvertex==nvertex4proc(myobj)) THEN
              myvertex = nvertex4proc(myobj-1)+1
           ENDIF

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

           ! if the vertex is the last one of the object, then it's also the first one
           IF(vertexindexlocal(mylocalvertex)==nvertex4proc(myobj)) THEN
              ! run another loop to find the 2nd vertex
              DO mylocalvertex1=nvertexlocal(myobj-1)+1,nvertexlocal(myobj)
                 IF(vertexindexlocal(mylocalvertex1)==nvertex4proc(myobj-1)+1) THEN

                    ! use the jc on 2nd vertex to minus the last(first) vertex
                    dudxdt = (ujc(1,mylocalvertex1+1) - ujc(1,mylocalvertex))/dist
                    dvdxdt = (vjc(1,mylocalvertex1+1) - vjc(1,mylocalvertex))/dist
                    dpdxdt = (pjc(1,mylocalvertex1+1) - pjc(1,mylocalvertex))/dist

                    dudydt = (ujc(2,mylocalvertex1+1) - ujc(2,mylocalvertex))/dist
                    dvdydt = (vjc(2,mylocalvertex1+1) - vjc(2,mylocalvertex))/dist
                    dpdydt = (pjc(2,mylocalvertex1+1) - pjc(2,mylocalvertex))/dist

                 ENDIF
              ENDDO
           ELSE
              ! if not the last vertex, just next minus current
              dudxdt = (ujc(1,mylocalvertex+1) - ujc(1,mylocalvertex))/dist
              dvdxdt = (vjc(1,mylocalvertex+1) - vjc(1,mylocalvertex))/dist
              dpdxdt = (pjc(1,mylocalvertex+1) - pjc(1,mylocalvertex))/dist

              dudydt = (ujc(2,mylocalvertex+1) - ujc(2,mylocalvertex))/dist
              dvdydt = (vjc(2,mylocalvertex+1) - vjc(2,mylocalvertex))/dist
              dpdydt = (pjc(2,mylocalvertex+1) - pjc(2,mylocalvertex))/dist

           ENDIF

           ! dduxjc dduyjc dduxyjc
           dduxjc=(tx*dudxdt-ty*dudydt+ty*ty*deltaujc1(mylocalvertex))/gacobi2
           ddvxjc=(tx*dvdxdt-ty*dvdydt+ty*ty*deltavjc1(mylocalvertex))/gacobi2

           dduyjc=(ty*dudydt-tx*dudxdt+tx*tx*deltaujc1(mylocalvertex))/gacobi2
           ddvyjc=(ty*dvdydt-tx*dvdxdt+tx*tx*deltavjc1(mylocalvertex))/gacobi2

           dduxyjc=(ty*dudxdt+tx*dudydt-tx*ty*deltaujc1(mylocalvertex))/gacobi2
           ddvxyjc=(ty*dvdxdt+tx*dvdydt-tx*ty*deltavjc1(mylocalvertex))/gacobi2

           ! ddpdx ddpdy jump conditions
           r3 = 0.0d0
           ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           r3 = 1.0d0
           ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2


           ujc(3,mylocalvertex)=dduxjc
           vjc(3,mylocalvertex)=ddvxjc
           pjc(3,mylocalvertex)=ddpxjc1
           ujc(5,mylocalvertex)=dduxyjc
           vjc(5,mylocalvertex)=ddvxyjc
           pjc(7,mylocalvertex)=ddpxjc0

           ujc(4,mylocalvertex)=dduyjc
           vjc(4,mylocalvertex)=ddvyjc
           pjc(4,mylocalvertex)=ddpyjc1
           ujc(6,mylocalvertex)=dduxyjc
           vjc(6,mylocalvertex)=ddvxyjc
           pjc(8,mylocalvertex)=ddpyjc0
        ENDIF
     ENDDO
  ENDDO

  ! clockwise direction
  DO myobj=1,nobj4proc
     DO mylocalvertex=nvertexlocal(myobj),nvertexlocal(myobj-1)+1,-1
        !get index of vertices
        myvertex=vertexindexlocal(mylocalvertex)
        IF(vertexindexlocal(mylocalvertex)==nvertex4proc(myobj-1)+1) THEN
           myvertex=nvertex4proc(myobj)
        ENDIF

        ! tao direction
        tx= -taox(myvertex-1)
        ty= -taoy(myvertex-1)

        gacobi=SQRT(tx*tx+ty*ty)
        gacobi2=tx*tx+ty*ty

        ! normal direction
        xn=-ty
        yn= tx
        xn=xn/gacobi
        yn=yn/gacobi

        ! velocity and position of vertices
        uu(0)=us(myvertex)
        vv(0)=vs(myvertex)
        xx(0)=vertex(1,myvertex)
        yy(0)=vertex(2,myvertex)
        distance(0)=zero

        ! dudn jump conditions on negative side
        ! dudnjcn= thetat X n

        dudnjcn = -thetat(myobj)*yn
        dvdnjcn =  thetat(myobj)*xn

        ! dudn jump conditions on positive side
        ! ============ put this part into a subroutine ============
        ! a. find 3 points on normal direction s1, s2, s3

        DO n=1,nany

           ! coordinates of interpolation point
           xx(n)= xi2x( x2xi( vertex(1,myvertex)) +DBLE(n)*ds*xn)
           yy(n)= eta2y(y2eta(vertex(2,myvertex)) +DBLE(n)*ds*yn)
           distance(n)=SQRT((xx(n)-xx(0))*(xx(n)-xx(0))+(yy(n)-yy(0))*(yy(n)-yy(0)))

           ! indices of interpolation point
           i=INT((x2xi(xx(n))-x2xi(x0))/dxi*2)
           j=INT((y2eta(yy(n))-y2eta(y0))/deta*2)
           IF(MOD(i,2)==0) THEN
              ie=i/2
              ic=ie
           ELSE
              ic=(i+1)/2
              ie=ic-1
           ENDIF
           IF(MOD(j,2)==0) THEN
              je=j/2
              jc=je
           ELSE
              jc=(j+1)/2
              je=jc-1
           ENDIF
           iu=ie
           ju=jc
           iv=ic
           jv=je

           id=INT(SIGN(1.d0,xn))
           jd=INT(SIGN(1.d0,yn))
           IF(id<0.0d0) THEN
              iu=iu+1
              iv=iv+1
           ENDIF
           IF(jd<0.0d0) THEN
              ju=ju+1
              jv=jv+1
           ENDIF

           ! interpolation of uu(n)
           DO j=1,many
              DO i=1,many
                 xa(i)=xf(iu+id*(i-1))
                 ya(i)=u(iu+id*(i-1),ju+jd*(j-1))
              ENDDO
              xb(j)=yc(ju+jd*(j-1))
              CALL interpolate(xa,ya,many,xx(n),yb(j),foo)
           ENDDO
           CALL interpolate(xb,yb,many,yy(n),uu(n),foo)

           ! interpolation of vv(n)
           DO i=1,many
              DO j=1,many
                 xa(j)=yf(jv+jd*(j-1))
                 ya(j)=v(iv+id*(i-1),jv+jd*(j-1))
              ENDDO
              xb(i)=xc(iv+id*(i-1))
              CALL interpolate(xa,ya,many,yy(n),yb(i),foo)
           ENDDO
           CALL interpolate(xb,yb,many,xx(n),vv(n),foo)
        ENDDO

        ! b. dudn jump conditions on positive side
        CALL FiniteDiff(distance,uu,dudnjcp,ddudnp)
        CALL FiniteDiff(distance,vv,dvdnjcp,ddvdnp)

        ! c. dudn jump conditions
        dudnjc = dudnjcp - dudnjcn
        dvdnjc = dvdnjcp - dvdnjcn

        ! laplacian u jump conditions
        ! ddudn jump conditions
        deltaujc2(mylocalvertex) = ddudnp + vertex(3,myvertex)*dudnjc
        deltavjc2(mylocalvertex) = ddvdnp + vertex(3,myvertex)*dvdnjc

        ! dudx dvdx dudy dvdy
        duxjc = -ty/(tx*yn-ty*xn)*dudnjc
        dvxjc = -ty/(tx*yn-ty*xn)*dvdnjc
        duyjc =  tx/(tx*yn-ty*xn)*dudnjc
        dvyjc =  tx/(tx*yn-ty*xn)*dvdnjc

        ! dpdx dpdy jump conditions
        ! body force qx, qy
        qx =  thetatt(myobj)*(vertex(2,myvertex)-ysc(myobj))
        qy = -thetatt(myobj)*(vertex(1,myvertex)-xsc(myobj))

        ! dpdx, dpdy jump conditions
        dpxjc = 1.0d0/Re*deltaujc1(mylocalvertex) + qx
        dpyjc = 1.0d0/re*deltavjc1(mylocalvertex) + qy

        ujc2(1,mylocalvertex)=duxjc
        vjc2(1,mylocalvertex)=dvxjc
        pjc2(1,mylocalvertex)=dpxjc
        ! print*, mylocalvertex,myvertex,ujc2(1,mylocalvertex),vjc2(1,mylocalvertex),pjc2(1,mylocalvertex)
        ujc2(2,mylocalvertex)=duyjc
        vjc2(2,mylocalvertex)=dvyjc
        pjc2(2,mylocalvertex)=dpyjc

     ENDDO

     ! second order jump conditions
     DO mylocalvertex=nvertexlocal(myobj),nvertexlocal(myobj-1)+1,-1
        myvertex=vertexindexlocal(mylocalvertex)
        A=vertex(:,myvertex)
        corner0=allcorner0(:,myid)
        corner1=allcorner1(:,myid)

        ! verify if this vertex in this domain
        IF( (A(1)>corner0(1)).AND.(A(2)>corner0(2)).AND. &
             (A(1)<corner1(1)).AND.(A(2)<corner1(2))) THEN
           IF(myvertex==nvertex4proc(myobj-1)+1) THEN
              myvertex = nvertex4proc(myobj)
           ENDIF

           ! tao direction
           tx=-taox(myvertex-1)
           ty=-taoy(myvertex-1)
           gacobi2=tx*tx+ty*ty

           ! dudxdt, dvdxdt, dudydt, dxdydt
           ! dudxdt = (dudxB-dudxA)/|AB|
           ! position of neighbor vertex
           xs0=vertex(1,myvertex)
           xs1=vertex(1,myvertex-1)
           ys0=vertex(2,myvertex)
           ys1=vertex(2,myvertex-1)
           dist=SQRT((ys1-ys0)*(ys1-ys0)+(xs1-xs0)*(xs1-xs0))

           ! if the vertex is the first one of the object, then it's also the last one
           IF(vertexindexlocal(mylocalvertex)==nvertex4proc(myobj-1)+1) THEN
              ! run another loop to find the 2nd last vertex
              DO mylocalvertex1=nvertexlocal(myobj-1)+1,nvertexlocal(myobj)
                 IF(vertexindexlocal(mylocalvertex1)==nvertex4proc(myobj)) THEN

                    ! use the jc on 2nd last vertex to minus the last(first) vertex
                    dudxdt = (ujc2(1,mylocalvertex1-1) - ujc2(1,mylocalvertex))/dist
                    dvdxdt = (vjc2(1,mylocalvertex1-1) - vjc2(1,mylocalvertex))/dist
                    dpdxdt = (pjc2(1,mylocalvertex1-1) - pjc2(1,mylocalvertex))/dist

                    dudydt = (ujc2(2,mylocalvertex1-1) - ujc2(2,mylocalvertex))/dist
                    dvdydt = (vjc2(2,mylocalvertex1-1) - vjc2(2,mylocalvertex))/dist
                    dpdydt = (pjc2(2,mylocalvertex1-1) - pjc2(2,mylocalvertex))/dist

                 ENDIF
              ENDDO
           ELSE
              ! if not the last vertex, just next minus current
              dudxdt = (ujc2(1,mylocalvertex-1) - ujc2(1,mylocalvertex))/dist
              dvdxdt = (vjc2(1,mylocalvertex-1) - vjc2(1,mylocalvertex))/dist
              dpdxdt = (pjc2(1,mylocalvertex-1) - pjc2(1,mylocalvertex))/dist

              dudydt = (ujc2(2,mylocalvertex-1) - ujc2(2,mylocalvertex))/dist
              dvdydt = (vjc2(2,mylocalvertex-1) - vjc2(2,mylocalvertex))/dist
              dpdydt = (pjc2(2,mylocalvertex-1) - pjc2(2,mylocalvertex))/dist

           ENDIF

           ! dduxjc dduyjc dduxyjc
           dduxjc=(tx*dudxdt-ty*dudydt+ty*ty*deltaujc2(mylocalvertex))/gacobi2
           ddvxjc=(tx*dvdxdt-ty*dvdydt+ty*ty*deltavjc2(mylocalvertex))/gacobi2

           dduyjc=(ty*dudydt-tx*dudxdt+tx*tx*deltaujc2(mylocalvertex))/gacobi2
           ddvyjc=(ty*dvdydt-tx*dvdxdt+tx*tx*deltavjc2(mylocalvertex))/gacobi2

           dduxyjc=(ty*dudxdt+tx*dudydt-tx*ty*deltaujc2(mylocalvertex))/gacobi2
           ddvxyjc=(ty*dvdxdt+tx*dvdydt-tx*ty*deltavjc2(mylocalvertex))/gacobi2

           ! ddpdx ddpdy jump conditions
           r3 = 0.0d0
           ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

           r3 = 1.0d0
           ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
           ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2


           ujc2(3,mylocalvertex)=dduxjc
           vjc2(3,mylocalvertex)=ddvxjc
           pjc2(3,mylocalvertex)=ddpxjc1
           ujc2(5,mylocalvertex)=dduxyjc
           vjc2(5,mylocalvertex)=ddvxyjc
           pjc2(7,mylocalvertex)=ddpxjc0

           ujc2(4,mylocalvertex)=dduyjc
           vjc2(4,mylocalvertex)=ddvyjc
           pjc2(4,mylocalvertex)=ddpyjc1
           ujc2(6,mylocalvertex)=dduxyjc
           vjc2(6,mylocalvertex)=ddvxyjc
           pjc2(8,mylocalvertex)=ddpyjc0
        ENDIF
     ENDDO
  ENDDO

  DO myobj=1,nobj4proc
     DO mylocalvertex=nvertexlocal(myobj-1)+1,nvertexlocal(myobj)

        ujc(1,mylocalvertex)=(ujc(1,mylocalvertex)+ujc2(1,mylocalvertex))/2.0d0
        ujc(2,mylocalvertex)=(ujc(2,mylocalvertex)+ujc2(2,mylocalvertex))/2.0d0
        ujc(3,mylocalvertex)=(ujc(3,mylocalvertex)+ujc2(3,mylocalvertex))/2.0d0
        ujc(4,mylocalvertex)=(ujc(4,mylocalvertex)+ujc2(4,mylocalvertex))/2.0d0
        ujc(5,mylocalvertex)=(ujc(5,mylocalvertex)+ujc2(5,mylocalvertex))/2.0d0
        ujc(6,mylocalvertex)=(ujc(6,mylocalvertex)+ujc2(6,mylocalvertex))/2.0d0

        vjc(1,mylocalvertex)=(vjc(1,mylocalvertex)+vjc2(1,mylocalvertex))/2.0d0
        vjc(2,mylocalvertex)=(vjc(2,mylocalvertex)+vjc2(2,mylocalvertex))/2.0d0
        vjc(3,mylocalvertex)=(vjc(3,mylocalvertex)+vjc2(3,mylocalvertex))/2.0d0
        vjc(4,mylocalvertex)=(vjc(4,mylocalvertex)+vjc2(4,mylocalvertex))/2.0d0
        vjc(5,mylocalvertex)=(vjc(5,mylocalvertex)+vjc2(5,mylocalvertex))/2.0d0
        vjc(6,mylocalvertex)=(vjc(6,mylocalvertex)+vjc2(6,mylocalvertex))/2.0d0

        pjc(1,mylocalvertex)=(pjc(1,mylocalvertex)+pjc2(1,mylocalvertex))/2.0d0
        pjc(2,mylocalvertex)=(pjc(2,mylocalvertex)+pjc2(2,mylocalvertex))/2.0d0
        pjc(3,mylocalvertex)=(pjc(3,mylocalvertex)+pjc2(3,mylocalvertex))/2.0d0
        pjc(4,mylocalvertex)=(pjc(4,mylocalvertex)+pjc2(4,mylocalvertex))/2.0d0
        pjc(5,mylocalvertex)=(pjc(5,mylocalvertex)+pjc2(5,mylocalvertex))/2.0d0
        pjc(6,mylocalvertex)=(pjc(6,mylocalvertex)+pjc2(6,mylocalvertex))/2.0d0
        pjc(7,mylocalvertex)=(pjc(7,mylocalvertex)+pjc2(7,mylocalvertex))/2.0d0
        pjc(8,mylocalvertex)=(pjc(8,mylocalvertex)+pjc2(8,mylocalvertex))/2.0d0

!if(myid == 2 .and. mylocalvertex<10) print*, mylocalvertex, vertexindexlocal(mylocalvertex), pjc(1,mylocalvertex)
     ENDDO
  ENDDO

  DEALLOCATE(deltaujc1)
  DEALLOCATE(deltavjc1)
  DEALLOCATE(deltaujc2)
  DEALLOCATE(deltavjc2)

  DEALLOCATE(ujc2)
  DEALLOCATE(vjc2)
  DEALLOCATE(pjc2)



END SUBROUTINE jc_firstsecond
