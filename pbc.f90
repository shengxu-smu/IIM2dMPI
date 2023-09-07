SUBROUTINE pbc(fac)

  USE Euler
  USE Lagrange
  USE para

  DOUBLE PRECISION:: conv,visc,dpdx,dpdy,fac

  !c  dirichlet
  IF(ipbcx0==1) THEN
     IF(ipbeg==1) THEN
        DO j=jpbeg,jpend
           rhsp(ipbeg,j)=COS(xc(i))*SIN(yc(j))
        ENDDO
     ENDIF
  ENDIF

  IF(ipbcx1==1) THEN
     IF(ipend==nx) THEN
        DO j=jpbeg,jpend
           rhsp(ipend,j)=COS(xc(i))*SIN(yc(j))
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy0==1) THEN
     IF(jpbeg==1) THEN
        DO i=ipbeg,ipend
           rhsp(i,jpbeg)=SIN(xc(i))*COS(yc(j))
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy1==1) THEN
     IF(jpend==ny) THEN
        DO i=ipbeg,ipend
           rhsp(i,jpend)=SIN(xc(i))*COS(yc(j))
        ENDDO
     ENDIF
  ENDIF

  !c  newmann
  IF(ipbcx0==2) THEN
     IF(ipbeg==1) THEN
        i=1
        DO j=jpbeg,jpend
           dpdx= COS(xc(i))*SIN(yc(j))
           rhsp(i,j)=rhsp(i,j)+2.0d0*xifx(i-1)*dxi*dpdx
        ENDDO
     ENDIF
  ENDIF

  IF(ipbcx1==2) THEN
     IF(ipend==nx) THEN
        i=nx
        DO j=jpbeg,jpend
           dpdx= COS(xc(i))*SIN(yc(j))
           rhsp(i,j)=rhsp(i,j)-2.0d0*xifx(i)*dxi*dpdx
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy0==2) THEN
     IF(jpbeg==1) THEN
        j=1
        DO i=ipbeg,ipend
           dpdy= SIN(xc(i))*COS(yc(j))
           rhsp(i,j)=rhsp(i,j)+2.0d0*etafy(j-1)*dxi*dxi/deta*dpdy
        ENDDO
     ENDIF
  ENDIF

  IF(jpbcy1==2) THEN
     IF(jpend==ny) THEN
        j=ny
        DO i=ipbeg,ipend
           dpdy= SIN(xc(i))*COS(yc(j))
           rhsp(i,j)=rhsp(i,j)-2.0d0*etafy(j)*dxi*dxi/deta*dpdy
        ENDDO
     ENDIF
  ENDIF

END SUBROUTINE pbc
