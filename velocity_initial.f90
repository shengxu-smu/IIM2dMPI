SUBROUTINE velocity_initial

  USE mpi
  USE para
  USE Euler
  USE Lagrange

  DOUBLE PRECISION::grad,ucpdx,vcpdy
  INTEGER:: ixc,jyc
  CALL exchange4ghost
  CALL interpolateuv(0)

  ! initialize rhsp
  grad = zero
  rhsp = zero
  IF(isingular.EQV..TRUE.) THEN
     ! find intersection points and allocate related jcs
     CALL surface_property
     CALL Raycrossing
     CALL allocate_jcs
     ! compute necessary jcs
    !  CALL jc_firstsecond
    !  CALL principle_jc_p
    !  CALL pressure_distribute

     ! compute jump contributions
    !  DO index = 1, n_xc_int
    !     CALL correction_pressure(index,1,vcpdy)
    !     ixc = xc_int_vertex(2, index)
    !     jyc = xc_int_vertex(3, index)
    !     v(ixc,jyc)=v(ixc,jyc)-vcpdy
    !  ENDDO
    !  DO index = 1, n_yc_int
    !     CALL correction_pressure(index,3,ucpdx)
    !     jyc = yc_int_vertex(2, index)
    !     ixc = yc_int_vertex(3, index)
    !     u(ixc,jyc)=u(ixc,jyc)-ucpdx
    !  ENDDO
  ENDIF

  ! solve pressure
!  CALL solvep

  ! update u & v
  DO j=jubeg,juend
    !  PRINT*,j,u(0,j)
     DO i=iubeg+1,iuend-1
        grad=-(p(i+1,j)-p(i,j))*xifx(i)/dxi
        u(i,j)=u(i,j)+grad
     ENDDO
  ENDDO
  DO j=jvbeg+1,jvend-1
     DO i=ivbeg,ivend
        grad=-(p(i,j+1)-p(i,j))*etafy(j)/deta
        v(i,j)=v(i,j)+grad
     ENDDO
  ENDDO

  IF(isingular.EQV..TRUE.) THEN
     CALL velocity_reset
  ENDIF

  CALL exchange4ghost
  CALL interpolateuv(0)
  IF(object_move.EQV..TRUE.) CALL free_jcs

END SUBROUTINE velocity_initial
