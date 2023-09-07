SUBROUTINE mac_distribute

  USE Euler
  USE para
  USE Lagrange

  DOUBLE PRECISION, DIMENSION(1:2)::A,B,corner0,corner1
  DOUBLE PRECISION:: signx,signy,f(22)
  DOUBLE PRECISION:: xsm0,xsm1,ysm0,ysm1,xx,yy
  INTEGER:: int_idx, myvertex

  ! xc
  DO int_idx = 1, n_xc_int
     myvertex = xc_int_vertex(1,int_idx)
     xx = xc_int_info(1, int_idx)
     !  yy = xc_int_info(2, int_idx)
     xsm0 = vertex(1,myvertex)
     xsm1 = vertex(1,myvertex+1)
     ysm0 = vertex(2,myvertex)
     ysm1 = vertex(2,myvertex+1)

     DO n=1,6
        ujcxc(n,int_idx)=ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
        vjcxc(n,int_idx)=vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
     ENDDO
    !  PRINT*,myvertex,vjc(2,myvertex),vjc(4,myvertex),vjc(6,myvertex)
    !  PRINT*,myvertex,vjcxc(2,int_idx),vjcxc(4,int_idx),vjcxc(6,int_idx)
     DO n=1,8
        pjcxc(n,int_idx)=pjc(n,myvertex)+(pjc(n,myvertex+1)-pjc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
     ENDDO
  ENDDO

  ! xf
  DO int_idx = 1,n_xf_int
     myvertex = xf_int_vertex(1,int_idx)
     xx = xf_int_info(1, int_idx)
     signx = xf_int_info(3,i)
     signy = xf_int_info(4,i)
     !  yy = xc_int_info(2, int_idx)
     xsm0 = vertex(1,myvertex)
     xsm1 = vertex(1,myvertex+1)
     ysm0 = vertex(2,myvertex)
     ysm1 = vertex(2,myvertex+1)
     DO n=1,6
        ujcxf(n,int_idx) = ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
        vjcxf(n,int_idx) = vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (xx-xsm0)/(xsm1-xsm0)
     ENDDO
     !print*,myvertex,ujcxf(1,int_idx),ujcxf(3,int_idx),ujcxf(5,int_idx)
  ENDDO

  ! yc
  DO int_idx = 1, n_yc_int
     myvertex = yc_int_vertex(1,int_idx)
     yy = yc_int_info(1, int_idx)
     !  xx = yc_int_info(2, int_idx)
     xsm0 = vertex(1,myvertex)
     xsm1 = vertex(1,myvertex+1)
     ysm0 = vertex(2,myvertex)
     ysm1 = vertex(2,myvertex+1)
     DO n=1,6
        ujcyc(n,int_idx)=ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
        vjcyc(n,int_idx)=vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
     ENDDO
     DO n=1,8
        pjcyc(n,int_idx)=pjc(n,myvertex)+(pjc(n,myvertex+1)-pjc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
     ENDDO
     !  if(myid == 0) PRINT*, int_idx, myvertex,xx,yy
    !  print*,myvertex,ujcyc(1,int_idx),ujcyc(3,int_idx),ujcyc(5,int_idx)
  ENDDO

  ! yf
  DO int_idx = 1,n_yf_int
     myvertex = yf_int_vertex(1,int_idx)
     yy = yf_int_info(1, int_idx)
     xsm0 = vertex(1,myvertex)
     xsm1 = vertex(1,myvertex+1)
     ysm0 = vertex(2,myvertex)
     ysm1 = vertex(2,myvertex+1)
     DO n=1,6
        ujcyf(n,int_idx)=ujc(n,myvertex)+(ujc(n,myvertex+1)-ujc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
        vjcyf(n,int_idx)=vjc(n,myvertex)+(vjc(n,myvertex+1)-vjc(n,myvertex))*&
             (yy-ysm0)/(ysm1-ysm0)
     ENDDO
  ENDDO

  !xc
  DO i = 1, n_xc_int
     signx = xc_int_info(3,i)
     signy = xc_int_info(4,i)

     myvertex=xc_int_vertex(1,i)
     !print*,myvertex,signx,signy
     ujcxc(1,i)=signx*ujcxc(1,i)
     vjcxc(1,i)=signx*vjcxc(1,i)
     ujcxc(3,i)=signx*ujcxc(3,i)
     vjcxc(3,i)=signx*vjcxc(3,i)
     ujcxc(5,i)=signx*ujcxc(5,i)
     vjcxc(5,i)=signx*vjcxc(5,i)

     ujcxc(2,i)=signy*ujcxc(2,i)
     vjcxc(2,i)=signy*vjcxc(2,i)
     ujcxc(4,i)=signy*ujcxc(4,i)
     vjcxc(4,i)=signy*vjcxc(4,i)
     ujcxc(6,i)=signy*ujcxc(6,i)
     vjcxc(6,i)=signy*vjcxc(6,i)
    !  PRINT*,myvertex,vjcxc(1,i),vjcxc(3,i),vjcxc(5,i)
  ENDDO

  ! xf
  DO i = 1, n_xf_int
     signx = xf_int_info(3,i)
     signy = xf_int_info(4,i)

     ujcxf(1,i)=signx*ujcxf(1,i)
     vjcxf(1,i)=signx*vjcxf(1,i)
     ujcxf(3,i)=signx*ujcxf(3,i)
     vjcxf(3,i)=signx*vjcxf(3,i)
     ujcxf(5,i)=signx*ujcxf(5,i)
     vjcxf(5,i)=signx*vjcxf(5,i)

     ujcxf(2,i)=signy*ujcxf(2,i)
     vjcxf(2,i)=signy*vjcxf(2,i)
     ujcxf(4,i)=signy*ujcxf(4,i)
     vjcxf(4,i)=signy*vjcxf(4,i)
     ujcxf(6,i)=signy*ujcxf(6,i)
     vjcxf(6,i)=signy*vjcxf(6,i)
  ENDDO

  ! yc
  DO j = 1, n_yc_int
     signx = yc_int_info(3,j)
     signy = yc_int_info(4,j)
     myvertex=yc_int_vertex(1,j)
     ujcyc(1,j)=signx*ujcyc(1,j)
     vjcyc(1,j)=signx*vjcyc(1,j)
     ujcyc(3,j)=signx*ujcyc(3,j)
     vjcyc(3,j)=signx*vjcyc(3,j)
     ujcyc(5,j)=signx*ujcyc(5,j)
     vjcyc(5,j)=signx*vjcyc(5,j)

     ujcyc(2,j)=signy*ujcyc(2,j)
     vjcyc(2,j)=signy*vjcyc(2,j)
     ujcyc(4,j)=signy*ujcyc(4,j)
     vjcyc(4,j)=signy*vjcyc(4,j)
     ujcyc(6,j)=signy*ujcyc(5,j)
     vjcyc(6,j)=signy*vjcyc(6,j)
     !print*,myvertex,ujcyc(1,j),ujcyc(3,j),ujcyc(5,j)
  ENDDO

  ! yf
  DO j = 1, n_yf_int
     signx = yf_int_info(3,j)
     signy = yf_int_info(4,j)

     ujcyf(1,j)=signx*ujcyf(1,j)
     vjcyf(1,j)=signx*vjcyf(1,j)
     ujcyf(3,j)=signx*ujcyf(3,j)
     vjcyf(3,j)=signx*vjcyf(3,j)
     ujcyf(5,j)=signx*ujcyf(5,j)
     vjcyf(5,j)=signx*vjcyf(5,j)

     ujcyf(2,j)=signy*ujcyf(2,j)
     vjcyf(2,j)=signy*vjcyf(2,j)
     ujcyf(4,j)=signy*ujcyf(4,j)
     vjcyf(4,j)=signy*vjcyf(4,j)
     ujcyf(6,j)=signy*ujcyf(6,j)
     vjcyf(6,j)=signy*vjcyf(6,j)
  ENDDO




END SUBROUTINE mac_distribute
