SUBROUTINE Euler_reset(krk)

  USE Euler
  USE para

  INTEGER,INTENT(in):: krk

  d = zero
  rhsp = zero
  ux = zero
  uy = zero
  vx = zero
  vy = zero
  rhsu(:,:,krk) = zero
  rhsv(:,:,krk) = zero
  IF(krk.EQ.1) THEN
     rhsu1 = zero
     rhsv1 = zero
  ELSEIF ( krk.EQ.2 ) THEN
     rhsu2 = zero
     rhsv2 = zero
  ELSEIF ( krk.EQ.3 ) THEN
     rhsu3 = zero
     rhsv3 = zero
  ELSE
     rhsu4 = zero
     rhsv4 = zero
  ENDIF

END SUBROUTINE Euler_reset
