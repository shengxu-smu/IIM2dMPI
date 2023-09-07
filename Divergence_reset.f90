SUBROUTINE Divergence_reset

  USE para
  USE Euler
  use Lagrange

  INTEGER:: ixc,ixc1,jyc,jyc1,index
  DOUBLE PRECISION:: signx,signy

  DO j=jdbeg,jdend
     DO i=idbeg,idend
        DO myobj=1,nobj4proc
           IF(iop(i,j).EQ.myobj) THEN
              d(i,j)=zero
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  DO index = 1, n_xc_int
     signy = xc_int_info(4, index)
     ixc = xc_int_vertex(2, index)
     jyc = xc_int_vertex(3, index)
     jyc1 = jyc+1
     IF(signy.GE.0.0d0) THEN
        d(ixc,jyc1)=2.0d0*d(ixc,jyc1+1)-d(ixc,jyc1+2)
     ELSE
        d(ixc,jyc)=2.0d0*d(ixc,jyc-1)-d(ixc,jyc-2)
     ENDIF
  ENDDO

  DO index = 1, n_yc_int
     signx = yc_int_info(3, index)
     jyc = yc_int_vertex(2, index)
     ixc = yc_int_vertex(3, index)
     ixc1 = ixc+1
     IF(signx.GE.0.0d0) THEN
        d(ixc1,jyc)=2.0d0*d(ixc1+1,jyc)-d(ixc1+2,jyc)
     ELSE
        d(ixc,jyc)=2.0d0*d(ixc-1,jyc)-d(ixc-2,jyc)
     ENDIF
  ENDDO

END SUBROUTINE Divergence_reset
