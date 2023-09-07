SUBROUTINE old_save

  USE para
  USE Euler

  DO j=jubeg,juend
     DO i=iubeg,iuend
        un(i,j)=u(i,j)
     ENDDO
  ENDDO
  DO j=jvbeg,jvend
     DO i=ivbeg,ivend
        vn(i,j)=v(i,j)
     ENDDO
  ENDDO

  DO j=jpbeg-1,jpend+1
     DO i=ipbeg-1,ipend+1
        dn(i,j)=ux(i,j)+vy(i,j)
     ENDDO
  ENDDO

END SUBROUTINE old_save
