SUBROUTINE correction_velocity(index,flag,krk,dpdx)

  USE stretchxy
  USE para
  USE Lagrange
  USE Euler

  INTEGER:: ixc,ixc1,ixf,ixf1,jyc,jyc1,jyf,jyf1
  INTEGER,INTENT(in):: index,flag,krk
  DOUBLE PRECISION,INTENT(in):: dpdx
  DOUBLE PRECISION:: ucudx,vcdxx(2),vcdx,ucdxx(2)
  DOUBLE PRECISION:: vcvdy,vcdyy(2),ucdy,ucdyy(2)
  DOUBLE PRECISION:: xint,yint,signx,signy
  DOUBLE PRECISION:: ui,vi,djc,ddjc
  DOUBLE PRECISION:: ujcdxi,ujcdxi2,ujcdeta,ujcdeta2
  DOUBLE PRECISION:: vjcdxi,vjcdxi2,vjcdeta,vjcdeta2
  DOUBLE PRECISION:: ucpdx,vcpdy

  ! xc
  IF(flag == 1) THEN

     yint = xc_int_info(2, index)
     signx = xc_int_info(3, index)
     signy = xc_int_info(4, index)

     ixc = xc_int_vertex(2, index)
     jyc = xc_int_vertex(3, index)
     jyf = xc_int_vertex(4, index)

     ui=xc_int_info(5,index)
     vi=xc_int_info(6,index)

     jyf1=jyf+1
     jyc1=jyc+1

     vjcdeta = vjcxc(2,index)/detady(y2eta(yint))
     vjcdeta2 = (vjcxc(4,index)-d2etady2(y2eta(yint))*vjcxc(2,index))/detady(y2eta(yint))/detady(y2eta(yint))

     djc=2.0d0*vi*vjcxc(2,index)/detady(y2eta(yint))
     ddjc=2.0d0*vi*(vjcxc(4,index)-d2etady2(y2eta(yint))*vjcxc(2,index))/detady(y2eta(yint))/detady(y2eta(yint)) &
          +2.0d0*(xc_int_info(10,index)*xc_int_info(10,index)- &
          xc_int_info(11,index)*xc_int_info(11,index))*detady(y2eta(yint))*detady(y2eta(yint))

     IF(jyc == jyf) THEN

        vcvdy=-1.0d0/deta*(djc*(y2eta(yc(jyc1))-y2eta(yint))+ &
             0.5d0*ddjc*(y2eta(yc(jyc1))-y2eta(yint))**2.0d0)

        vcdyy(1)=(-detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (vjcdeta*(y2eta(yc(jyc1))-y2eta(yint))+0.5d0*vjcdeta2*(y2eta(yc(jyc1))-y2eta(yint))**2.0d0)
        vcdyy(2)=( detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (vjcdeta*(y2eta(yc(jyc)) -y2eta(yint))+0.5d0*vjcdeta2*(y2eta(yc(jyc)) -y2eta(yint))**2.0d0)
     ELSE

        vcvdy=-1.0d0/deta*(djc*(y2eta(yc(jyc))-y2eta(yint))+ &
             0.5d0*ddjc*(y2eta(yc(jyc))-y2eta(yint))**2.0d0)

        vcdyy(1)=(-detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (vjcdeta*(y2eta(yf(jyf1))-y2eta(yint))+0.5d0*vjcdeta2*(y2eta(yf(jyf1))-y2eta(yint))**2.0d0)
        vcdyy(2)=( detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (vjcdeta*(y2eta(yf(jyf)) -y2eta(yint))+0.5d0*vjcdeta2*(y2eta(yf(jyf)) -y2eta(yint))**2.0d0)

     ENDIF
     vcpdy = dpdx
     !  print*,ixc,jyc,rhsv1(ixc,jyc),rhsv1(ixc,jyf),rhsv1(ixc,jyf1)
    !  PRINT*,ixc,jyc,rhsv(ixc,jyc,krk),rhsv(ixc,jyf,krk),rhsv(ixc,jyf1,krk)
     rhsv(ixc,jyc, krk) = rhsv(ixc,jyc,krk)  - vcvdy - vcpdy
     rhsv(ixc,jyf, krk) = rhsv(ixc,jyf,krk)  + vcdyy(1)/Re
     rhsv(ixc,jyf1,krk) = rhsv(ixc,jyf1,krk) + vcdyy(2)/Re
     IF( krk.EQ.1 ) THEN
        rhsv1(ixc,jyc)  = rhsv1(ixc,jyc)  - vcvdy - vcpdy
        rhsv1(ixc,jyf)  = rhsv1(ixc,jyf)  + vcdyy(1)/Re
        rhsv1(ixc,jyf1) = rhsv1(ixc,jyf1) + vcdyy(2)/Re
     ELSEIF ( krk.EQ.2 ) THEN
        rhsv2(ixc,jyc)  = rhsv2(ixc,jyc)  - vcvdy - vcpdy
        rhsv2(ixc,jyf)  = rhsv2(ixc,jyf)  + vcdyy(1)/Re
        rhsv2(ixc,jyf1) = rhsv2(ixc,jyf1) + vcdyy(2)/Re
     ELSEIF ( krk.EQ.3 ) THEN
        rhsv3(ixc,jyc)  = rhsv3(ixc,jyc)  - vcvdy - vcpdy
        rhsv3(ixc,jyf)  = rhsv3(ixc,jyf)  + vcdyy(1)/Re
        rhsv3(ixc,jyf1) = rhsv3(ixc,jyf1) + vcdyy(2)/Re
     ELSE
        rhsv3(ixc,jyc)  = rhsv4(ixc,jyc)  - vcvdy - vcpdy
        rhsv3(ixc,jyf)  = rhsv4(ixc,jyf)  + vcdyy(1)/Re
        rhsv3(ixc,jyf1) = rhsv4(ixc,jyf1) + vcdyy(2)/Re
     ENDIF
    !  PRINT*,ixc,jyc,rhsv(ixc,jyc,krk),rhsv(ixc,jyf,krk),rhsv(ixc,jyf1,krk)
     ! print*,ixc,jyc,rhsv1(ixc,jyc),rhsv1(ixc,jyf),rhsv1(ixc,jyf1)
     ! print*,ixc,jyc,vcvdy,vcpdy

  ENDIF

  ! xf
  IF(flag == 2) THEN

     yint = xf_int_info(2, index)
     signx = xf_int_info(3, index)
     signy = xf_int_info(4, index)

     ixf = xf_int_vertex(2, index)
     jyc = xf_int_vertex(3, index)
     jyf = xf_int_vertex(4, index)

     ui=xf_int_info(5,index)
     vi=xf_int_info(6,index)

     jyf1=jyf+1
     jyc1=jyc+1

     ujcdeta = ujcxf(2,index)/detady(y2eta(yint))
     ujcdeta2 = (ujcxf(4,index)-d2etady2(y2eta(yint))*ujcxf(2,index))/detady(y2eta(yint))/detady(y2eta(yint))

     djc=(ui*vjcxf(2,index)+vi*ujcxf(2,index))/detady(y2eta(yint))
     ddjc=ui*(vjcxf(4,index)-d2etady2(y2eta(yint))*vjcxf(2,index))/detady(y2eta(yint))/detady(y2eta(yint))+&
          vi*(ujcxf(4,index)-d2etady2(y2eta(yint))*ujcxf(2,index))/detady(y2eta(yint))/detady(y2eta(yint))+&
          2.0d0*(xf_int_info(7,index)*xf_int_info(9,index)- &
          xf_int_info(8,index)*xf_int_info(10,index))*detady(y2eta(yint))*detady(y2eta(yint))

     IF(jyc == jyf) THEN
        ucdy=-1.0d0/deta*(djc*(y2eta(yf(jyf))-y2eta(yint))+ &
             0.5d0*ddjc*(y2eta(yf(jyf))-y2eta(yint))**2.0d0)

        ucdyy(1)=(-detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (ujcdeta*(y2eta(yc(jyc1))-y2eta(yint))+0.5d0*ujcdeta2*(y2eta(yc(jyc1))-y2eta(yint))**2.0d0)
        ucdyy(2)=( detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (ujcdeta*(y2eta(yc(jyc)) -y2eta(yint))+0.5d0*ujcdeta2*(y2eta(yc(jyc)) -y2eta(yint))**2.0d0)

     ELSE
        ucdy=-1.0d0/deta*(djc*(y2eta(yf(jyf1))-y2eta(yint))+ &
             0.5d0*ddjc*(y2eta(yf(jyf1))-y2eta(yint))**2.0d0)

        ucdyy(1)=(-detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (ujcdeta*(y2eta(yf(jyf1))-y2eta(yint))+0.5d0*ujcdeta2*(y2eta(yf(jyf1))-y2eta(yint))**2.0d0)
        ucdyy(2)=( detady(y2eta(yint))*detady(y2eta(yint))/deta/deta-d2etady2(y2eta(yint))/2.0d0/deta)*&
             (ujcdeta*(y2eta(yf(jyf)) -y2eta(yint))+0.5d0*ujcdeta2*(y2eta(yf(jyf)) -y2eta(yint))**2.0d0)

     ENDIF
    !  print*,ixf,jyc,djc,ddjc
     IF( krk.EQ.1 ) THEN
        rhsu1(ixf,jyf1) = rhsu1(ixf,jyf1) - ucdy
        rhsu1(ixf,jyf)  = rhsu1(ixf,jyf)  + ucdyy(1)*Re1
        rhsu1(ixf,jyf1) = rhsu1(ixf,jyf1) + ucdyy(2)*Re1
     ELSEIF ( krk.EQ.2 ) THEN
        rhsu2(ixf,jyf1) = rhsu2(ixf,jyf1) - ucdy
        rhsu2(ixf,jyf)  = rhsu2(ixf,jyf)  + ucdyy(1)*Re1
        rhsu2(ixf,jyf1) = rhsu2(ixf,jyf1) + ucdyy(2)*Re1
     ELSEIF ( krk.EQ.3 ) THEN
        rhsu3(ixf,jyf1) = rhsu3(ixf,jyf1) - ucdy
        rhsu3(ixf,jyf)  = rhsu3(ixf,jyf)  + ucdyy(1)*Re1
        rhsu3(ixf,jyf1) = rhsu3(ixf,jyf1) + ucdyy(2)*Re1
     ELSE
        rhsu4(ixf,jyf1) = rhsu4(ixf,jyf1) - ucdy
        rhsu4(ixf,jyf)  = rhsu4(ixf,jyf)  + ucdyy(1)*Re1
        rhsu4(ixf,jyf1) = rhsu4(ixf,jyf1) + ucdyy(2)*Re1
     ENDIF
     rhsu(ixf,jyf1,krk) = rhsu(ixf,jyf1,krk) - ucdy
     rhsu(ixf,jyf,krk) = rhsu(ixf,jyf,krk) + ucdyy(1)*Re1
     rhsu(ixf,jyf1,krk) = rhsu(ixf,jyf1,krk) + ucdyy(2)*Re1

  ENDIF

  ! yc
  IF(flag == 3) THEN

     xint = yc_int_info(2, index)
     signx = yc_int_info(3, index)
     signy = yc_int_info(4, index)

     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)

     ui = yc_int_info(5,index)
     vi = yc_int_info(6,index)

     ixf1=ixf+1
     ixc1=ixc+1

     djc=2.0d0*ui*ujcyc(1,index)/dxidx(x2xi(xint))
     ddjc=2.0d0*ui*(ujcyc(3,index)-d2xidx2(x2xi(xint))*ujcyc(1,index))/dxidx(x2xi(xint))/dxidx(x2xi(xint))+&
          2.0d0*(yc_int_info(7,index)*yc_int_info(7,index)-yc_int_info(8,index)*yc_int_info(8,index)) &
          /dxidx(x2xi(xint))/dxidx(x2xi(xint))

     ujcdxi = ujcyc(1,index)/dxidx(x2xi(xint))
     ujcdxi2 = (ujcyc(3,index)-d2xidx2(x2xi(xint))*ujcyc(1,index))/dxidx(x2xi(xint))/dxidx(x2xi(xint))

     IF(ixc == ixf) THEN

        ucudx=-1.0d0/dxi*(djc*(x2xi(xc(ixc1))-x2xi(xint))+ &
             0.5d0*ddjc*(x2xi(xc(ixc1))-x2xi(xint))**2.0d0)
        ucdxx(1)=(-dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (ujcdxi*(x2xi(xf(ixf1))-x2xi(xint))+0.5d0*ujcdxi2*(x2xi(xf(ixf1))-x2xi(xint))**2.0d0)
        ucdxx(2)=( dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (ujcdxi*(x2xi(xf(ixf)) -x2xi(xint))+0.5d0*ujcdxi2*(x2xi(xf(ixf)) -x2xi(xint))**2.0d0)
     ELSE
        ucudx=-1.0d0/dxi*(djc*(x2xi(xc(ixc))-x2xi(xint))+ &
             0.5d0*ddjc*(x2xi(xc(ixc))-x2xi(xint))**2.0d0)
        ucdxx(1)=(-dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (ujcdxi*(x2xi(xf(ixf1))-x2xi(xint))+0.5d0*ujcdxi2*(x2xi(xf(ixf1))-x2xi(xint))**2.0d0)
        ucdxx(2)=( dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (ujcdxi*(x2xi(xf(ixf)) -x2xi(xint))+0.5d0*ujcdxi2*(x2xi(xf(ixf)) -x2xi(xint))**2.0d0)
     ENDIF
     ucpdx = dpdx
     IF( krk.EQ.1 ) THEN
        rhsu1(ixc,jyc) = rhsu1(ixc,jyc) - ucudx - ucpdx
        rhsu1(ixf,jyc) = rhsu1(ixf,jyc) + ucdxx(1)*Re1
        rhsu1(ixf1,jyc) = rhsu1(ixf1,jyc) + ucdxx(2)*Re1
     ELSEIF ( krk.EQ.2 ) THEN
        rhsu2(ixc,jyc) = rhsu2(ixc,jyc) - ucudx - ucpdx
        rhsu2(ixf,jyc) = rhsu2(ixf,jyc) + ucdxx(1)*Re1
        rhsu2(ixf1,jyc) = rhsu2(ixf1,jyc) + ucdxx(2)*Re1
     ELSEIF ( krk.EQ.3 ) THEN
        rhsu3(ixc,jyc) = rhsu3(ixc,jyc) - ucudx - ucpdx
        rhsu3(ixf,jyc) = rhsu3(ixf,jyc) + ucdxx(1)*Re1
        rhsu3(ixf1,jyc) = rhsu3(ixf1,jyc) + ucdxx(2)*Re1
     ELSE
        rhsu4(ixc,jyc) = rhsu4(ixc,jyc) - ucudx - ucpdx
        rhsu4(ixf,jyc) = rhsu4(ixf,jyc) + ucdxx(1)*Re1
        rhsu4(ixf1,jyc) = rhsu4(ixf1,jyc) + ucdxx(2)*Re1
     ENDIF
     rhsu(ixc,jyc,krk) = rhsu(ixc,jyc,krk) - ucudx - ucpdx
     rhsu(ixf,jyc,krk) = rhsu(ixf,jyc,krk) + ucdxx(1)*Re1
     rhsu(ixf1,jyc,krk) = rhsu(ixf1,jyc,krk) + ucdxx(2)*Re1

  ENDIF

  ! yf
  IF(flag == 4) THEN

     xint = yf_int_info(2, index)
     signx = yf_int_info(3, index)
     signy = yf_int_info(4, index)

     jyf = yf_int_vertex(2, index)
     ixc = yf_int_vertex(3, index)
     ixf = yf_int_vertex(4, index)

     ui = yf_int_info(5,index)
     vi = yf_int_info(6,index)

     ixf1=ixf+1
     ixc1=ixc+1

     djc=(ui*vjcyf(1,index)+vi*ujcyf(1,index))/dxidx(x2xi(xint))
     ddjc=ui*(vjcyf(3,index)-d2xidx2(x2xi(xint))*vjcyf(1,index))/dxidx(x2xi(xint))/dxidx(x2xi(xint))+&
          vi*(ujcyf(3,index)-d2xidx2(x2xi(xint))*ujcyf(1,index))/dxidx(x2xi(xint))/dxidx(x2xi(xint))+&
          2.0d0*(yf_int_info(7,index)*yf_int_info(9,index)- &
          yf_int_info(8,index)*yf_int_info(10,index))*dxidx(x2xi(xint))*x2xi(x2xi(xint))

     vjcdxi  = vjcyf(2,index)/dxidx(x2xi(xint))
     vjcdxi2 = (vjcyf(4,index)-d2xidx2(x2xi(xint))*vjcyf(2,index))/dxidx(x2xi(xint))/dxidx(x2xi(xint))

     IF(ixc == ixf) THEN
        vcdx=-1.0d0/dxi*(djc*(x2xi(xf(ixf))-x2xi(xint))+0.5d0*ddjc*(x2xi(xf(ixf))-x2xi(xint))**2.0d0)
        vcdxx(1)=(-dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (vjcdxi*(x2xi(xc(ixc1))-x2xi(xint))+0.5d0*vjcdxi2*(x2xi(xc(ixc1))-x2xi(xint))**2.0d0)
        vcdxx(2)=( dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (vjcdxi*(x2xi(xc(ixc)) -x2xi(xint))+0.5d0*vjcdxi2*(x2xi(xc(ixc)) -x2xi(xint))**2.0d0)

     ELSE
        vcdx=-1.0d0/dxi*(djc*(x2xi(xf(ixf1))-x2xi(xint))+0.5d0*ddjc*(x2xi(xf(ixf1))-x2xi(xint))**2.0d0)
        vcdxx(1)=(-dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (vjcdxi*(x2xi(xc(ixc1))-x2xi(xint))+0.5d0*vjcdxi2*(x2xi(xc(ixc1))-x2xi(xint))**2.0d0)
        vcdxx(2)=( dxidx(x2xi(xint))*dxidx(x2xi(xint))/dxi/dxi-d2xidx2(x2xi(xint))/2.0d0/dxi)*&
             (vjcdxi*(x2xi(xc(ixc)) -x2xi(xint))+0.5d0*vjcdxi2*(x2xi(xc(ixc)) -x2xi(xint))**2.0d0)
     ENDIF

     IF( krk.EQ.1 ) THEN
        rhsv1(ixf1,jyf) = rhsv1(ixf1,jyf) - vcdx
        rhsv1(ixc,jyf) = rhsv1(ixf,jyf) + vcdxx(1)*Re1
        rhsv1(ixc1,jyf) = rhsv1(ixc1,jyf) + vcdxx(2)*Re1
     ELSEIF ( krk.EQ.2 ) THEN
        rhsv1(ixf1,jyf) = rhsv1(ixf1,jyf) - vcdx
        rhsv1(ixc,jyf) = rhsv1(ixf,jyf) + vcdxx(1)*Re1
        rhsv1(ixc1,jyf) = rhsv1(ixc1,jyf) + vcdxx(2)*Re1
     ELSEIF ( krk.EQ.3 ) THEN
        rhsv1(ixf1,jyf) = rhsv1(ixf1,jyf) - vcdx
        rhsv1(ixc,jyf) = rhsv1(ixf,jyf) + vcdxx(1)*Re1
        rhsv1(ixc1,jyf) = rhsv1(ixc1,jyf) + vcdxx(2)*Re1
     ELSE
        rhsv1(ixf1,jyf) = rhsv1(ixf1,jyf) - vcdx
        rhsv1(ixc,jyf) = rhsv1(ixf,jyf) + vcdxx(1)*Re1
        rhsv1(ixc1,jyf) = rhsv1(ixc1,jyf) + vcdxx(2)*Re1
     ENDIF
     rhsv(ixf1,jyf,krk) = rhsv(ixf1,jyf,krk) - vcdx
     rhsv(ixc,jyf,krk) = rhsv(ixf,jyf,krk) + vcdxx(1)*Re1
     rhsv(ixc1,jyf,krk) = rhsv(ixc1,jyf,krk) + vcdxx(2)*Re1

  ENDIF


END SUBROUTINE correction_velocity
