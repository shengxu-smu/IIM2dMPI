SUBROUTINE rhs4p(fac)

  USE para
  USE Euler
  USE stretchxy

  DOUBLE PRECISION:: fac

  ! IF(isingular.EQV..TRUE.) CALL Divergence_reset
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        IF(iop(i,j).NE.0) THEN
           rhsp(i,j) = rhsp(i,j)+2.0d0*EXP(-xc(i))*EXP(-yc(j))*dxi*dxi
        ELSE
           rhsp(i,j) = rhsp(i,j)-2.0d0*SIN( xc(i))*SIN( yc(j))*dxi*dxi
        ENDIF
     ENDDO
  ENDDO
  CALL pbc

END SUBROUTINE rhs4p
