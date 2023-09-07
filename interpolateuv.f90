SUBROUTINE interpolateuv(flag)

  USE para
  USE Euler
  USE Lagrange

  INTEGER::flag
  DOUBLE PRECISION:: w1,w2

  IF(flag==0 .OR. flag==4) THEN
     ! uee vee
     DO j=jpbeg-1,jpend
        DO i=ipbeg-1,ipend
           w1=(yf(j)-yc(j))  /(yc(j+1)-yc(j))
           w2=(yc(j+1)-yf(j))/(yc(j+1)-yc(j))
           uee(i,j)= u(i,j)*w1+u(i,j+1)*w2

           w1=(xf(i)-xc(i))  /(xc(i+1)-xc(i))
           w2=(xc(i+1)-xf(i))/(xc(i+1)-xc(i))
           vee(i,j)= v(i,j)*w1+v(i+1,j)*w2
        ENDDO
     ENDDO
  ENDIF

  IF (flag==0 .OR. flag==3) THEN
     ! ucc
     DO j=jpbeg-1,jpend+1
        DO i=iubeg,iuend
           w1=(xc(i)-xf(i-1))/(xf(i)-xf(i-1))
           w2=(xf(i)-xc(i))  /(xf(i)-xf(i-1))
           ucc(i,j)= u(i-1,j)*w1+u(i,j)*w2
        ENDDO
     ENDDO
  ENDIF

  IF (flag==0 .OR. flag==2) THEN
     ! vec
     DO j=jpbeg,jpend
        DO i=ipbeg-1,ipend
           w1=(yc(j)-yf(j-1))/(yf(j)-yf(j-1))
           w2=(yf(j)-  yc(j))/(yf(j)-yf(j-1))
           vec(i,j)=vee(i,j-1)*w1+vee(i,j)*w2
        ENDDO
     ENDDO
  ENDIF

  IF (flag == 0 .OR. flag==1) THEN
     ! uce
     DO j=jpbeg-1,jpend
        DO i=ipbeg,ipend
           w1=(yf(j)-  yc(j))/(yc(j+1)-yc(j))
           w2=(yc(j+1)-yf(j))/(yc(j+1)-yc(j))
           uce(i,j)=ucc(i,j)*w1+ucc(i,j+1)*w2
        ENDDO
     ENDDO
     ! vcc
     DO j=jubeg,jvend
        DO i=ipbeg-1,ipend+1
           w1=(yc(j)-yf(j-1))/(yf(j)-yf(j-1))
           w2=(yf(j)-yc(j))  /(yf(j)-yf(j-1))
           vcc(i,j)= v(i,j-1)*w1+v(i,j)*w2
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE interpolateuv
