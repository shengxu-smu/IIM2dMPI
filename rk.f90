SUBROUTINE rk

  USE mpi
  USE Euler
  USE Lagrange
  USE para

  DOUBLE PRECISION:: fac
  INTEGER:: krk

  ! krk=1:
  krk=1
  fac=0.5d0
  CALL Euler_reset(krk)
  IF(isingular.EQV..TRUE.) THEN
     IF(object_move.EQV..TRUE.) THEN
        CALL surface_property
        CALL Raycrossing
        CALL allocate_jcs
     ENDIF
     CALL jc_firstsecond
     CALL principle_jc_p
     CALL mac_distribute
     CALL jump_contribution(krk)
  ELSE
    !  CALL interpolateuv(0)
  ENDIF
  ! CALL uv_strain
  ! CALL old_save
  CALL rhs4p(fac)
  CALL solvep
  ! CALL rhs4uv(krk)
  ! CALL updateuv(krk,fac)
  ! CALL exchange4ghost
  ! IF(object_move.EQV..TRUE.) CALL free_jcs

  ! krk=2:
  ! krk=2
  ! fac=0.5d0
  ! CALL Euler_reset(krk)
  ! IF(isingular.EQV..TRUE.) THEN
  !    IF(object_move.EQV..TRUE.) THEN
  !       CALL surface_property
  !       CALL Raycrossing
  !       CALL allocate_jcs
  !    ENDIF
  !    CALL jc_firstsecond
  !    CALL principle_jc_p
  !    CALL mac_distribute
  !    CALL jump_contribution(krk)
  ! ELSE
  !    CALL interpolateuv(0)
  ! ENDIF
  ! CALL uv_strain
  ! CALL rhs4p(fac)
  ! CALL solvep
  ! CALL rhs4uv(krk)
  ! CALL updateuv(krk,fac)
  ! CALL exchange4ghost
  ! IF(object_move.EQV..TRUE.) CALL free_jcs
  !
  ! ! krk=3:
  ! krk=3
  ! fac=1.0d0
  ! CALL Euler_reset(krk)
  ! IF(isingular.EQV..TRUE.) THEN
  !    IF(object_move.EQV..TRUE.) THEN
  !       CALL surface_property
  !       CALL Raycrossing
  !       CALL allocate_jcs
  !    ENDIF
  !    CALL jc_firstsecond
  !    CALL principle_jc_p
  !    CALL mac_distribute
  !    CALL jump_contribution(krk)
  ! ELSE
  !    CALL interpolateuv(0)
  ! ENDIF
  ! CALL uv_strain
  ! CALL rhs4p(fac)
  ! CALL solvep
  ! CALL rhs4uv(krk)
  ! CALL updateuv(krk,fac)
  ! CALL exchange4ghost
  ! IF(object_move.EQV..TRUE.) CALL free_jcs
  !
  ! ! krk=4:
  ! krk=4
  ! fac=1.0d0
  ! CALL Euler_reset(krk)
  ! IF(isingular.EQV..TRUE.) THEN
  !    IF(object_move.EQV..TRUE.) THEN
  !       CALL surface_property
  !       CALL Raycrossing
  !       CALL allocate_jcs
  !    ENDIF
  !    CALL jc_firstsecond
  !    CALL principle_jc_p
  !    CALL mac_distribute
  !    CALL jump_contribution(krk)
  ! ELSE
  !    CALL interpolateuv(0)
  ! ENDIF
  ! CALL uv_strain
  ! CALL rhs4p(fac)
  ! CALL solvep
  ! CALL rhs4uv(krk)
  ! CALL updateuv(krk,fac)
  ! CALL exchange4ghost
  ! IF(object_move.EQV..TRUE.) CALL free_jcs

END SUBROUTINE rk
