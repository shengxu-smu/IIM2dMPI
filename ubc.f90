SUBROUTINE ubc

  USE para
  USE euler

  DOUBLE PRECISION::w1,w2,ubcx0,ubcx1,ubcy0,ubcy1
  ! Periodic BCs


  ! Dirichlet BCs
  ubcx0=1.0d0
  ubcx1=0.0d0
  ubcy0=0.0d0
  ubcy1=0.0d0

  IF(iubcx0==1) THEN
     IF(iubeg==0) THEN
        DO j=jubeg,juend
           w1=(xc(1)-xf(0))/(xf(1)-xf(0))
           w2=(xf(1)-xc(1))/(xf(1)-xf(0))
           u(0,j)=(ubcx0-u(1,j)*w2)/w1
        ENDDO
     ENDIF
  ENDIF

  IF(iubcx1==1) THEN
     IF(iuend==nx) THEN
        DO j=jubeg,juend
           w1=(xc(nx)-xf(nx-1))/(xf(nx)-xf(nx-1))
           w2=(xf(nx)-xc(nx))  /(xf(nx)-xf(nx-1))
           u(nx,j)=(ubcx1-u(nx-1,j)*w1)/w2
        ENDDO
     ENDIF
  ENDIF

  IF(jubcy0==1) THEN
     IF(jubeg==0) THEN
        DO i=iubeg,iuend
           w1=(yc(1)-yc(0))/(yc(2)-yc(0))
           w2=(yc(2)-yc(1))/(yc(2)-yc(0))
           u(i,0)=(ubcy0-u(i,2)*w2)/w1
        ENDDO
     ENDIF
  ENDIF

  IF(jubcy1==1) THEN
     IF(juend==ny+1) THEN
        DO i=iubeg,iuend
           w1=(yc(ny)-yc(ny-1))/(yc(ny+1)-yc(ny-1))
           w2=(yc(ny+1)-yc(ny))/(yc(ny+1)-yc(ny-1))
           u(i,ny+1)=(1.0d0-u(i,ny-1)*w1)/w2
        ENDDO
     ENDIF
  ENDIF


  ! Neumann BCs
  IF(iubcx0==2) THEN
     IF(iubeg==0) THEN
        DO j=jubeg,juend
           u(0,j)=u(1,j)-0.0d0*(dxi/xicx(1))
        ENDDO
     ENDIF
  ENDIF

  IF(iubcx1==2) THEN
     IF(iuend==nx) THEN
        DO j=jubeg,juend
           u(nx,j)=u(nx-1,j)+0.0d0*(dxi/xicx(nx))
        ENDDO
     ENDIF
  ENDIF

  IF(jubcy0==2) THEN
     IF(jubeg==0) THEN
        DO i=iubeg,iuend
           u(i,0)=u(i,2)-0.0d0*(deta*2.0d0/etacy(1))
        ENDDO
     ENDIF
  ENDIF

  IF(jubcy1==2) THEN
     IF(juend==ny+1) THEN
        DO i=iubeg,iuend
           u(i,ny+1)=u(i,ny-1)+0.0d0*(2.0d0*deta/etacy(ny))
        ENDDO
     ENDIF
  ENDIF

END SUBROUTINE ubc
