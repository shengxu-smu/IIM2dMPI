SUBROUTINE correction_strain(index, flag)

  USE para
  USE Euler
  USE Lagrange
  USE stretchxy

  INTEGER:: flag,index,ixc,ixf,ixf1,jyc,jyf,jyf1
  DOUBLE PRECISION:: xint,yint,pcudy,pcvdy,pcudx,pcvdx

  ! xc
  IF(flag == 1) THEN
     yint = xc_int_info(2,index)
     ixc = xc_int_vertex(2,index)
     jyc = xc_int_vertex(3,index)
     jyf = xc_int_vertex(4,index)
     jyf1 = jyf+1

     IF(jyc==jyf) THEN
        !      pcudy(i,l)=-dy1*(ujcxc(2,i,l)*(yf(jyf)-yint)+
        ! .               0.5d0*ujcxc(4,i,l)*(yf(jyf)-yint)**2.0d0)
        !      pcvdy(i,l)=-dy1*(vjcxc(2,i,l)*(yf(jyf)-yint)+
        ! .               0.5d0*vjcxc(4,i,l)*(yf(jyf)-yint)**2.0d0)
        pcudy=-1.0d0/deta*(ujcxc(2,index)*(y2eta(yf(jyf))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(ujcxc(4,index)-d2etady2(y2eta(yint))*ujcxc(2,index)) * &
             (y2eta(yf(jyf))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
        pcvdy=-1.0d0/deta*(vjcxc(2,index)*(y2eta(yf(jyf))-y2eta(yint)/detady(y2eta(yint)))+ &
             0.5d0*(vjcxc(4,index)-d2etady2(y2eta(yint))*vjcxc(2,index)) * &
             (y2eta(yf(jyf))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
     ELSE
        !      pcudy(i,l)=-dy1*(ujcxc(2,i,l)*(yf(jyf1)-yint)+
        ! .               0.5d0*ujcxc(4,i,l)*(yf(jyf1)-yint)**2.0d0)
        !      pcvdy(i,l)=-dy1*(vjcxc(2,i,l)*(yf(jyf1)-yint)+
        ! .               0.5d0*vjcxc(4,i,l)*(yf(jyf1)-yint)**2.0d0)
        pcudy=-1.0d0/deta*(ujcxc(2,index)*(y2eta(yf(jyf1))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(ujcxc(4,index)-d2etady2(y2eta(yint))*ujcxc(2,index)) * &
             (y2eta(yf(jyf1))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
        pcvdy=-1.0d0/deta*(vjcxc(2,index)*(y2eta(yf(jyf1))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(vjcxc(4,index)-d2etady2(y2eta(yint))*vjcxc(2,index)) * &
             (y2eta(yf(jyf1))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
     ENDIF
     uy(ixc,jyf1) = uy(ixc,jyf1) + pcudy
     vy(ixc,jyf1) = vy(ixc,jyf1) + pcvdy
     d(ixc,jyf1)  = d(ixc,jyf1)  + pcvdy
  ENDIF

  ! yc
  IF(flag == 3) THEN
     xint = yc_int_info(2, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)
     ixf1=ixf+1
     IF(ixc==ixf) THEN
        !      pcudx(j,l)=-dx1*(ujcyc(1,j,l)*(xf(ixf)-xint)+
        ! .               0.5d0*ujcyc(3,j,l)*(xf(ixf)-xint)**2.0d0)
        !      pcvdx(j,l)=-dx1*(vjcyc(1,j,l)*(xf(ixf)-xint)+
        ! .               0.5d0*vjcyc(3,j,l)*(xf(ixf)-xint)**2.0d0)
        pcudx=-1.0d0/dxi*(ujcyc(1,index)*(x2xi(xf(ixf))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(ujcyc(3,index)-d2xidx2(x2xi(xint))*ujcyc(1,index)) * &
             (x2xi(xf(ixf))-x2xi(xint))**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
        pcvdx=-1.0d0/dxi*(vjcyc(1,index)*(x2xi(xf(ixf))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(vjcyc(3,index)-d2xidx2(x2xi(xint))*vjcyc(1,index)) * &
             (x2xi(xf(ixf))-x2xi(xint))**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
     ELSE
        !      pcudx(j,l)=-dx1*(ujcyc(1,j,l)*(xf(ixf1)-xint)+
        ! .               0.5d0*ujcyc(3,j,l)*(xf(ixf1)-xint)**2.0d0)
        !      pcvdx(j,l)=-dx1*(vjcyc(1,j,l)*(xf(ixf1)-xint)+
        ! .               0.5d0*vjcyc(3,j,l)*(xf(ixf1)-xint)**2.0d0)
        pcudx=-1.0d0/dxi*(ujcyc(1,index)*(x2xi(xf(ixf1))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(ujcyc(3,index)-d2xidx2(x2xi(xint))*ujcyc(1,index)) * &
             (x2xi(xf(ixf1))-x2xi(xint))**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
        pcvdx=-1.0d0/dxi*(vjcyc(1,index)*(x2xi(xf(ixf1))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(vjcyc(3,index)-d2xidx2(x2xi(xint))*vjcyc(1,index)) * &
             (x2xi(xf(ixf1))-x2xi(xint))**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
     ENDIF
    !  print*,ixf1,jyc,ux(ixf1,jyc)
     ux(ixf1, jyc) = ux(ixf1,jyc) + pcudx
    !  print*,ixf1,jyc,ux(ixf1,jyc),pcudx
     vx(ixf1, jyc) = vx(ixf1,jyc) + pcvdx
     d(ixf1, jyc)  = d(ixf1, jyc) + pcudx

  ENDIF


END SUBROUTINE correction_strain
