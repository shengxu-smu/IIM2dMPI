SUBROUTINE rk

  USE mpi
  USE Euler
  USE Lagrange
  USE para

  DOUBLE PRECISION,DIMENSION(1:2)::a,b,c,cof1,cof2,cof3,cof0,dd
  DOUBLE PRECISION:: conv,visc,dpdx,dpdy,fac,summ1,summ2,ddd
  INTEGER:: krk
  ! RK3 for du/dt = rhsu
  ! where, rhsu = -A(u)-G(p)+D(u)
  !        A(u) = A1(u)+A2(u)+A3(u)
  !        D(u) = L(u)/Re

  ! krk=1:

  krk=1
  fac=0.5d0

  CALL exchange4ghost
  CALL interpolateuv
  CALL old_save
  CALL rhs4p(fac)
  CALL solvep
  CALL rhs4uv(krk)
  CALL updateuv(krk,fac)

  ! krk=2:
  ! krk=2
  ! fac=0.5d0
  !
  ! CALL exchange4ghost
  ! CALL interpolateuv
  ! CALL rhs4p(fac)
  ! CALL solvep
  ! CALL rhs4uv(krk)
  ! CALL updateuv(krk,fac)
  !
  ! ! krk=3:
  ! krk=3
  ! fac=1.0d0
  !
  ! CALL exchange4ghost
  ! CALL interpolateuv
  ! CALL rhs4p(fac)
  ! CALL solvep
  ! CALL rhs4uv(krk)
  ! CALL updateuv(krk,fac)
  !
  ! krk=4
  ! fac=1.0d0
  !
  ! CALL exchange4ghost
  ! CALL interpolateuv
  ! CALL rhs4p(fac)
  ! CALL solvep
  ! CALL rhs4uv(krk)
  ! CALL updateuv(krk,fac)

  CALL exchange4ghost
  CALL interpolateuv



END SUBROUTINE rk
