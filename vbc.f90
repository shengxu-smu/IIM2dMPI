SUBROUTINE vbc

  USE para
  USE Euler

  DOUBLE PRECISION::w1,w2,vbcx0,vbcx1,vbcy0,vbcy1
  ! dirichlet

  vbcx0=0.0d0
  vbcx1=0.0d0
  vbcy0=0.0d0
  vbcy1=0.0d0

  IF(ivbcx0==1) THEN
     IF(ivbeg==0) THEN
        DO j=jvbeg,jvend
           w1=(xc(1)-xc(0))/(xc(2)-xc(0))
           w2=(xc(2)-xc(1))/(xc(2)-xc(0))
           v(0,j)=(vbcx0-v(2,j)*w2)/w1
        ENDDO
     ENDIF
  ENDIF

  IF(ivbcx1==1) THEN
     IF(ivend==nx+1) THEN
        DO j=jvbeg,jvend
           w1=(xc(nx)-xc(nx-1))/(xc(nx+1)-xc(nx-1))
           w2=(xc(nx+1)-xc(nx))/(xc(nx+1)-xc(nx-1))
           v(nx+1,j)=(vbcx1-v(nx-1,j)*w1)/w2
        ENDDO
     ENDIF
  ENDIF

  IF(jvbcy0==1) THEN
     IF(jvbeg==0) THEN
        DO i=ivbeg,ivend
           w1=(yc(1)-yf(0))/(yf(1)-yf(0))
           w2=(yf(1)-yc(1))/(yf(1)-yf(0))
           v(i,0)=(vbcy0-v(i,1)*w2)/w1
        ENDDO
     ENDIF
  ENDIF

  IF(jvbcy1==1) THEN
     IF(jvend==ny) THEN
        DO i=ivbeg,ivend
           w1=(yc(ny)-yf(ny-1))/(yf(ny)-yf(ny-1))
           w2=(yf(ny)-yc(ny  ))/(yf(ny)-yf(ny-1))
           v(i,ny)=(vbcy1-v(i,ny-1)*w1)/w2
        ENDDO
     ENDIF
  ENDIF

  ! neumann
  IF(ivbcx0==2) THEN
     IF(ivbeg==0) THEN
        DO j=jvbeg,jvend
           v(0,j)=v(2,j)-0.0d0*(2.0d0*dxi/xicx(1))
        ENDDO
     ENDIF
  ENDIF

  IF(ivbcx1==2) THEN
     IF(ivend==nx+1) THEN
        DO j=jvbeg,jvend
           v(nx+1,j)=v(nx-1,j)+0.0d0*(2.0d0*dxi/xicx(nx))
        ENDDO
     ENDIF
  ENDIF

  IF(jvbcy0==2) THEN
     IF(jvbeg==0) THEN
        DO i=ivbeg,ivend
           v(i,0)=v(i,1)-0.0d0*(deta/etacy(1))
        ENDDO
     ENDIF
  ENDIF

  IF(jvbcy1==2) THEN
     IF(jvend==ny) THEN
        DO i=ivbeg,ivend
           v(i,ny)=v(i,ny-1)+0.0d0*(deta/etacy(ny))
        ENDDO
     ENDIF
  ENDIF


END SUBROUTINE vbc
