SUBROUTINE velocity_reset

  USE para
  USE Euler
  USE Lagrange

  DOUBLE PRECISION xx,yy

  DO j=jubeg,juend
     DO i=iubeg,iuend
        DO myobj=1,nobj4proc
           IF(iou(i,j).EQ.myobj) THEN
              u(i,j)=zero
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  DO j=jvbeg,jvend
     DO i=ivbeg,ivend
        DO myobj=1,nobj4proc
           IF(iov(i,j).EQ.myobj) THEN
              v(i,j)=zero
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE velocity_reset
