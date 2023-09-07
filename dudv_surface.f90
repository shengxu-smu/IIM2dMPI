SUBROUTINE dudv_surface(index, flag)

  USE stretchxy
  USE Euler
  USE para
  USE Lagrange

  INTEGER:: index,flag,ixc,ixf,jyc,jyf,myvertex
  DOUBLE PRECISION:: signx,signy
  DOUBLE PRECISION:: duxjc,duyjc,dvxjc,dvyjc
  DOUBLE PRECISION:: duyp,duyn,dvyp,dvyn,duxp,duxn,dvxp,dvxn

  DOUBLE PRECISION:: gacobi2,gacobi,tx,ty,xn,yn
  DOUBLE PRECISION:: dudnjcp,dvdnjcp,dump

  ! variables used in interpolation of dudnjc ddudnjc
  DOUBLE PRECISION:: ds,dist,foo
  DOUBLE PRECISION,DIMENSION(0:3):: uu,vv,xx,yy,distance
  INTEGER:: ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
  INTEGER, PARAMETER:: many=3,nany=3
  DOUBLE PRECISION,DIMENSION(1:many):: xa,ya,xb,yb

  ds=1.01d0*SQRT(dxi*dxi+deta*deta)
  ! xc
  IF(flag == 1) THEN

     signx = xc_int_info(3, index)
     signy = xc_int_info(4, index)

     myvertex = xc_int_vertex(1, index)
     ixc = xc_int_vertex(2, index)
     jyc = xc_int_vertex(3, index)
     jyf = xc_int_vertex(4, index)

     ! jump conditions
     duxjc=signy*signx*ujcxc(1,index)
     dvxjc=signy*signx*vjcxc(1,index)

     duyjc=ujcxc(2,index)
     dvyjc=vjcxc(2,index)

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
     uu(0)=xc_int_info(5, index)
     vv(0)=xc_int_info(6, index)

     xx(0)=xc_int_info(1, index)
     yy(0)=xc_int_info(2, index)
     distance(0)=zero

     ! dudn jump conditions on positive side
     ! ============ put this part into a subroutine ============
     ! a. find 3 points on normal direction s1, s2, s3

     DO n=1,nany

        ! coordinates of interpolation point
        xx(n)= xi2x( x2xi( xx(0)) +DBLE(n)*ds*xn)
        yy(n)= eta2y(y2eta(yy(0)) +DBLE(n)*ds*yn)
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
     CALL FiniteDiff(distance,uu,dudnjcp,dump)
     CALL FiniteDiff(distance,vv,dvdnjcp,dump)

     ! dudx dvdx dudy dvdy
     duxp = -ty/(tx*yn-ty*xn)*dudnjcp
     dvxp = -ty/(tx*yn-ty*xn)*dvdnjcp
     duyp =  tx/(tx*yn-ty*xn)*dudnjcp
     dvyp =  tx/(tx*yn-ty*xn)*dvdnjcp

     duxn = duxp - duxjc
     dvxn = dvxp - dvxjc
     duyn = duyp - duyjc
     dvyn = dvyp - dvyjc

     xc_int_info(7, index)=duxp
     xc_int_info(8, index)=duxn
     xc_int_info(9, index)=dvxp
     xc_int_info(10,index)=dvxn
     xc_int_info(11,index)=duyp
     xc_int_info(12,index)=duyn
     xc_int_info(13,index)=dvyp
     xc_int_info(14,index)=dvyn

  ENDIF

  ! xf
  IF(flag == 2) THEN

     signx = xf_int_info(3, index)
     signy = xf_int_info(4, index)

     myvertex = xf_int_vertex(1, index)
     ixf = xf_int_vertex(2, index)
     jyc = xf_int_vertex(3, index)
     jyf = xf_int_vertex(4, index)
     myobj = xf_int_vertex(5, index)

     ! jump conditions
     duyjc = ujcxf(2,index)
     dvyjc = vjcxf(2,index)

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
     uu(0)=xf_int_info(5, index)
     vv(0)=xf_int_info(6, index)

     xx(0)=xf_int_info(1, index)
     yy(0)=xf_int_info(2, index)
     distance(0)=zero

     ! dudn jump conditions on positive side
     ! ============ put this part into a subroutine ============
     ! a. find 3 points on normal direction s1, s2, s3

     DO n=1,nany

        ! coordinates of interpolation point
        xx(n)= xi2x( x2xi( xx(0)) +DBLE(n)*ds*xn)
        yy(n)= eta2y(y2eta(yy(0)) +DBLE(n)*ds*yn)
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
     CALL FiniteDiff(distance,uu,dudnjcp,dump)
     CALL FiniteDiff(distance,vv,dvdnjcp,dump)

     ! dudx dvdx dudy dvdy
     duyp = tx/(tx*yn-ty*xn)*dudnjcp
     dvyp = tx/(tx*yn-ty*xn)*dvdnjcp

     duyn = duyp - duyjc
     dvyn = dvyp - dvyjc

     xf_int_info(7, index)=duyp
     xf_int_info(8, index)=duyn
     xf_int_info(9, index)=dvyp
     xf_int_info(10,index)=dvyn

  ENDIF

  ! yc
  IF(flag == 3) THEN

     signx = yc_int_info(3, index)
     signy = yc_int_info(4, index)

     myvertex = yc_int_vertex(1, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)
     myobj = yc_int_vertex(5, index)

     ! jump conditions
     duxjc=ujcyc(1,index)
     dvxjc=vjcyc(1,index)

     duyjc=signx*signy*ujcyc(2,index)
     dvyjc=signx*signy*vjcyc(2,index)

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
     uu(0)=yc_int_info(5, index)
     vv(0)=yc_int_info(6, index)

     xx(0)=yc_int_info(2, index)
     yy(0)=yc_int_info(1, index)
     distance(0)=zero

     ! dudn jump conditions on positive side
     ! ============ put this part into a subroutine ============
     ! a. find 3 points on normal direction s1, s2, s3

     DO n=1,nany

        ! coordinates of interpolation point
        xx(n)= xi2x( x2xi( xx(0)) +DBLE(n)*ds*xn)
        yy(n)= eta2y(y2eta(yy(0)) +DBLE(n)*ds*yn)
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
     CALL FiniteDiff(distance,uu,dudnjcp,dump)
     CALL FiniteDiff(distance,vv,dvdnjcp,dump)

     ! dudx dvdx dudy dvdy
     duxp = -ty/(tx*yn-ty*xn)*dudnjcp
     dvxp = -ty/(tx*yn-ty*xn)*dvdnjcp
     duyp =  tx/(tx*yn-ty*xn)*dudnjcp
     dvyp =  tx/(tx*yn-ty*xn)*dvdnjcp

     duxn = duxp - duxjc
     dvxn = dvxp - dvxjc
     duyn = duyp - duyjc
     dvyn = dvyp - dvyjc

     yc_int_info(7, index)=duxp
     yc_int_info(8, index)=duxn
     yc_int_info(9, index)=dvxp
     yc_int_info(10,index)=dvxn
     yc_int_info(11,index)=duyp
     yc_int_info(12,index)=duyn
     yc_int_info(13,index)=dvyp
     yc_int_info(14,index)=dvyn

  ENDIF

  ! yf
  IF(flag == 4) THEN

     signx = yf_int_info(3, index)
     signy = yf_int_info(4, index)

     myvertex = yf_int_vertex(1, index)
     jyf = yf_int_vertex(2, index)
     ixc = yf_int_vertex(3, index)
     ixf = yf_int_vertex(4, index)
     myobj = yf_int_vertex(5, index)

     ! jump conditions
     duxjc = ujcyf(1,index)
     dvxjc = vjcyf(1,index)

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
     uu(0)=yf_int_info(5, index)
     vv(0)=yf_int_info(6, index)

     xx(0)=yf_int_info(2, index)
     yy(0)=yf_int_info(1, index)
     distance(0)=zero

     ! dudn jump conditions on positive side
     ! ============ put this part into a subroutine ============
     ! a. find 3 points on normal direction s1, s2, s3

     DO n=1,nany

        ! coordinates of interpolation point
        xx(n)= xi2x( x2xi( xx(0)) +DBLE(n)*ds*xn)
        yy(n)= eta2y(y2eta(yy(0)) +DBLE(n)*ds*yn)
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
     CALL FiniteDiff(distance,uu,dudnjcp,dump)
     CALL FiniteDiff(distance,vv,dvdnjcp,dump)

     ! dudx dvdx dudy dvdy
     duxp = -ty/(tx*yn-ty*xn)*dudnjcp
     dvxp = -ty/(tx*yn-ty*xn)*dvdnjcp

     duxn = duxp - duxjc
     dvxn = dvxp - dvxjc

     yf_int_info(7, index)=duxp
     yf_int_info(8, index)=duxn
     yf_int_info(9, index)=dvxp
     yf_int_info(10,index)=dvxn
     !  print*,ixf,jyf,yf_int_info(7,index),yf_int_info(8,index)
     ! print*,ixf,jyf,duxp,duxn,duxjc
  ENDIF

END SUBROUTINE dudv_surface
