SUBROUTINE correction_pressure(index,flag,dpdx)

  USE para
  USE Euler
  USE Lagrange
  USE stretchxy

  INTEGER:: ixf,ixc,ixc1,jyf,jyc,jyc1,myvertex
  INTEGER, INTENT(in):: flag,index
  DOUBLE PRECISION:: xint,yint,pcdyy(2),vcpdy,pcdxx(2),ucpdx
  DOUBLE PRECISION:: pjc0,pjcdxi,pjcdxi2,pjcdeta,pjcdeta2
  DOUBLE PRECISION, INTENT(out):: dpdx

  dpdx = zero
  IF(flag == 1) THEN
     myvertex=xc_int_vertex(1,index)
     ixc = xc_int_vertex(2,index)
     jyc = xc_int_vertex(3,index)
     jyf = xc_int_vertex(4,index)
     yint = xc_int_info(2,index)
     jyc1=jyc+1

     ! jump conditions of pressure are found in non-uniform grid, so need to be changed to uniform grid
     ! [p(xi)]=[p(x)]
     pjc0 = pjcxc(6,index)
     ! [dpdxi]=[dpdx]*dxdxi
     pjcdeta = pjcxc(2,index)*dydeta(yint)
     ! [d2pdxi2] = [d2pdx2]*dxdxi*dxdxi+[dpdx]*d2xdxi2
     pjcdeta2 = pjcxc(4,index)*dydeta(yint)*dydeta(yint) + pjcxc(2,index)*d2ydeta2(yint)

     ! in order to find jump contributions of pressure to non-uniform grid points, first we need to find
     ! them in uniform grids, then
     ! [d2pdx2] = dxidx*dxidx*[d2pdxi2] + d2xidxx*[dpdxi], where [d2pdxi2] and [dpdxi] are calculated from
     ! the formula on 2006_SISC_Xu_Wang.pdf
     pcdyy(1)=(-detady(yint)*detady(yint)/deta/deta - d2etady2(yint)/2.0d0/deta)*&
          (pjc0 + pjcdeta*(etac(jyc1)-y2eta(yint))  + 0.5d0*pjcdeta2*(etac(jyc1)-y2eta(yint))**2.0d0)
     pcdyy(2)=( detady(yint)*detady(yint)/deta/deta - d2etady2(yint)/2.0d0/deta)*&
          (pjc0 + pjcdeta*(etac(jyc)- y2eta(yint))  + 0.5d0*pjcdeta2*(etac(jyc)- y2eta(yint))**2.0d0)

     IF(jyc.EQ.jyf) THEN
        vcpdy=-1.0d0/deta/dydeta(yint)*(pjc0+ &
             pjcdeta *(etac(jyc1)-y2eta(yint))+ &
             0.5d0*pjcdeta2*(etac(jyc1)-y2eta(yint))**2.0d0)
     ELSE
        vcpdy=-1.0d0/deta/dydeta(yint)*(pjc0+ &
             pjcdeta *(etac(jyc)- y2eta(yint))+ &
             0.5d0*pjcdeta2*(etac(jyc)- y2eta(yint))**2.0d0)
     ENDIF
     ! update rhsp
     dpdx = vcpdy
     IF(ixc<=ipend .AND. ixc>=ipbeg) THEN
        IF(jyc <=jpend .AND. jyc >=jpbeg) rhsp(ixc,jyc)  =rhsp(ixc,jyc)  -pcdyy(1)*dxi*dxi
        IF(jyc1<=jpend .AND. jyc1>=jpbeg) rhsp(ixc,jyc+1)=rhsp(ixc,jyc+1)-pcdyy(2)*dxi*dxi
     ENDIF
  ENDIF

  IF(flag == 3) THEN
     myvertex=yc_int_vertex(1,index)
     xint = yc_int_info(2, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)
     ixc1=ixc+1

     pjc0 = pjcyc(5,index)
     pjcdxi = pjcyc(1,index)*dxdxi(xint)
     pjcdxi2 = pjcyc(3,index)*dxdxi(xint)*dxdxi(xint) + pjcyc(1,index)*d2xdxi2(xint)

     pcdxx(1)=(-dxidx(xint)*dxidx(xint)/dxi/dxi - d2xidx2(xint)/2.0d0/dxi)*&
          (pjc0 + pjcdxi*(xic(ixc1)-x2xi(xint)) + 0.5d0*pjcdxi2*(xic(ixc1)-x2xi(xint))**2.0d0)
     pcdxx(2)=( dxidx(xint)*dxidx(xint)/dxi/dxi - d2xidx2(xint)/2.0d0/dxi)*&
          (pjc0 + pjcdxi*(xic(ixc)- x2xi(xint)) + 0.5d0*pjcdxi2*(xic(ixc)- x2xi(xint))**2.0d0)

     IF(ixc.EQ.ixf) THEN
        ucpdx=-1.0d0/dxi/dxdxi(xint)*(pjc0+ &
             pjcdxi *(xic(ixc1)-x2xi(xint))+ &
             0.5d0*pjcdxi2*(xic(ixc1)-x2xi(xint))**2.0d0)
     ELSE
        ucpdx=-1.0d0/dxi/dxdxi(xint)*(pjc0+ &
             pjcdxi *(xic(ixc)- x2xi(xint))+ &
             0.5d0*pjcdxi2*(xic(ixc)- x2xi(xint))**2.0d0)
     ENDIF
     ! update rhsp
     IF(jyc>=jpbeg .AND. jyc<=jpend) THEN
        IF(ixc <=ipend .AND. ixc >=ipbeg) rhsp(ixc,  jyc)=rhsp(ixc,  jyc)-pcdxx(1)*dxi*dxi
        IF(ixc1<=ipend .AND. ixc1>=ipbeg) rhsp(ixc+1,jyc)=rhsp(ixc+1,jyc)-pcdxx(2)*dxi*dxi
     ENDIF
     dpdx = ucpdx
  ENDIF

END SUBROUTINE correction_pressure
