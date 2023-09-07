SUBROUTINE correction_interpolate(index, flag)

  USE stretchxy
  USE para
  USE Euler
  USE Lagrange

  ! indices
  INTEGER, INTENT(in):: flag,index
  INTEGER:: is,js,myvertex,ixc,ixf,jyc,jyf
  INTEGER:: u_jc_index_y, u_jc_index_x, v_jc_index_y,v_jc_index_x
  ! jump contribution
  DOUBLE PRECISION:: vcixfyf,vciycxc,vciycxf
  DOUBLE PRECISION:: uciyfxc,uciyfxf,ucixcyc
  DOUBLE PRECISION:: xint,yint

  ! xc
  IF(flag == 1) THEN
     ixc = xc_int_vertex(2,index)
     jyc = xc_int_vertex(3,index)
     jyf = xc_int_vertex(4,index)
     yint = xc_int_info(2,index)
     IF( jyc == jyf ) THEN
        v_jc_index_y = jyc+1
        js = jyf
        vciycxc = 0.5d0*(vjcxc(2,index)*(y2eta(yf(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(vjcxc(4,index)-d2etady2(y2eta(yint))*vjcxc(2,index)) * &
             (y2eta(yf(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
        u_jc_index_y = jyf
        js = jyc+1
        uciyfxc =-0.5d0*(ujcxc(2,index)*(y2eta(yc(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(ujcxc(4,index)-d2etady2(y2eta(yint))*ujcxc(2,index)) * &
             (y2eta(yc(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
     ELSE
        v_jc_index_y = jyc
        js = jyf+1
        vciycxc =-0.5d0*(vjcxc(2,index)*(y2eta(yf(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(vjcxc(4,index)-d2etady2(y2eta(yint))*vjcxc(2,index)) * &
             (y2eta(yf(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
        u_jc_index_y = jyf+1
        js = jyc
        uciyfxc = 0.5d0*(ujcxc(2,index)*(y2eta(yc(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(ujcxc(4,index)-d2etady2(y2eta(yint))*ujcxc(2,index)) * &
             (y2eta(yc(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
     ENDIF
     !  PRINT*,ixc, u_jc_index_y,uce(ixc, u_jc_index_y),ucc(ixc, u_jc_index_y),ucc(ixc, u_jc_index_y+1)
     vcc(ixc, v_jc_index_y) = vcc(ixc, v_jc_index_y) + vciycxc
     uce(ixc, u_jc_index_y) = uce(ixc, u_jc_index_y) + uciyfxc

  ENDIF
  ! xf
  IF(flag == 2) THEN
     ixf = xf_int_vertex(2,index)
     jyc = xf_int_vertex(3,index)
     jyf = xf_int_vertex(4,index)
     yint = xf_int_info(2,index)
     IF( jyc == jyf) THEN
        u_jc_index_y = jyf
        js = jyc+1
        ! uciyfxf = -0.5d0*(ujcxf(2,index)*(y2eta(yc(js))-y2eta(yint))/detady(y2eta(yint))+ &
        !      0.5d0*ujcxf(4,index)*(yc(js)-yint)**2.0d0)
        uciyfxf =-0.5d0*(ujcxf(2,index)*(y2eta(yc(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(ujcxf(4,index)-d2etady2(y2eta(yint))*ujcxf(2,index)) * &
             (y2eta(yc(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
        v_jc_index_y = jyc+1
        js = jyf
        vciycxf = 0.5d0*(vjcxf(2,index)*(y2eta(yf(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(vjcxf(4,index)-d2etady2(y2eta(yint))*vjcxf(2,index)) * &
             (y2eta(yf(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
        !  vciycxf = 0.5d0*(vjcxf(2,index)*(yf(js)-yint)+ &
        ! 0.5d0*vjcxf(4,index)*(yf(js)-yint)**2.0d0)
     ELSE
        u_jc_index_y = jyf+1
        js = jyc
        ! uciyfxf =0.5d0*(ujcxf(2,index)*(yc(js)-yint)+ &
        !      0.5d0*ujcxf(4,index)*(yc(js)-yint)**2.0d0)
        uciyfxf = 0.5d0*(ujcxf(2,index)*(y2eta(yc(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(ujcxf(4,index)-d2etady2(y2eta(yint))*ujcxf(2,index)) * &
             (y2eta(yc(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
        v_jc_index_y =jyc
        js = jyf+1
        ! vciycxf =-0.5d0*(vjcxf(2,index)*(yf(js)-yint)+ &
        !      0.5d0*vjcxf(4,index)*(yf(js)-yint)**2.0d0)
        vciycxf =-0.5d0*(vjcxf(2,index)*(y2eta(yf(js))-y2eta(yint))/detady(y2eta(yint))+ &
             0.5d0*(vjcxf(4,index)-d2etady2(y2eta(yint))*vjcxf(2,index)) * &
             (y2eta(yf(js))-y2eta(yint))**2.0d0/detady(y2eta(yint))/detady(y2eta(yint)))
     ENDIF
     !  print*,ixf,v_jc_index_y,vee(ixf, v_jc_index_y)
     uee(ixf, u_jc_index_y) = uee(ixf, u_jc_index_y) + uciyfxf
     vec(ixf, v_jc_index_y) = vec(ixf, v_jc_index_y) + vciycxf

  ENDIF
  ! yc
  IF(flag == 3) THEN
     xint = yc_int_info(2, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixf = yc_int_vertex(4, index)
     IF(ixc == ixf) THEN
        u_jc_index_x = ixc+1
        is = ixf
        ! ucixcyc=0.5d0*(ujcyc(1,index)*(xf(is)-xint)+ &
        !  0.5d0*ujcyc(3,index)*(xf(is)-xint)**2.0d0)
        ucixcyc= 0.5d0*(ujcyc(1,index)*(x2xi(xf(is))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(ujcyc(3,index)-d2xidx2(x2xi(xint))*ujcyc(1,index)) * &
             (x2xi(xf(is))-x2xi(xint))**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
     ELSE
        u_jc_index_x = ixc
        is = ixf+1
        ! ucixcyc=-0.5d0*(ujcyc(1,index)*(xf(is)-xint)+ &
        !  0.5d0*ujcyc(3,index)*(xf(is)-xint)**2.0d0)
        ucixcyc=-0.5d0*(ujcyc(1,index)*(x2xi(xf(is))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(ujcyc(3,index)-d2xidx2(x2xi(xint))*ujcyc(1,index)) * &
             (x2xi(xf(is))-x2xi(xint))**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
     ENDIF
     !  print*,u_jc_index_x,jyc,ucc(u_jc_index_x,jyc),ucixcyc
     ucc(u_jc_index_x, jyc) = ucc(u_jc_index_x, jyc) + ucixcyc
     !  print*,u_jc_index_x,jyc,ucc(u_jc_index_x,jyc),ucixcyc
  ENDIF

  ! yf
  IF(flag == 4) THEN
     xint = yf_int_info(2, index)
     jyf = yf_int_vertex(2, index)
     ixc = yf_int_vertex(3, index)
     ixf = yf_int_vertex(4, index)
     IF( ixc == ixf ) THEN
        v_jc_index_x = ixf
        is = ixc+1
        ! vcixfyf=-0.5d0*(vjcyf(1,index)*(xc(is)-xint)+ &
        !      0.5d0*vjcyf(3,index)*(xc(is)-xint)**2.0d0)
        vcixfyf=-0.5d0*(vjcyf(1,index)*(x2xi(xc(is))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(vjcyf(3,index)-d2xidx2(x2xi(xint))*vjcyf(1,index)) * &
             (xc(is)-xint)**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
     ELSE
        v_jc_index_x = ixf+1
        is = ixc
        ! vcixfyf= 0.5d0*(vjcyf(1,index)*(xc(is)-xint)+ &
        !      0.5d0*vjcyf(3,index)*(xc(is)-xint)**2.0d0)
        vcixfyf= 0.5d0*(vjcyf(1,index)*(x2xi(xc(is))-x2xi(xint))/dxidx(x2xi(xint))+ &
             0.5d0*(vjcyf(3,index)-d2xidx2(x2xi(xint))*vjcyf(1,index)) * &
             (xc(is)-xint)**2.0d0/dxidx(x2xi(xint))/dxidx(x2xi(xint)))
     ENDIF
     vee(v_jc_index_x,jyf) = vee(v_jc_index_x,jyf)+vcixfyf
  ENDIF

END SUBROUTINE correction_interpolate
