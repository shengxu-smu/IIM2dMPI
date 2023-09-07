SUBROUTINE uv_strain

  USE stretchxy
  USE para
  USE Euler
  use Lagrange

  ! DO index = 1, n_yc_int
  !    jyc = yc_int_vertex(2, index)
  !    ixc = yc_int_vertex(3, index)
  !    ixf = yc_int_vertex(4, index)
  !    ixf1=ixf+1
  !    print*,ixf1,jyc,ux(ixf1,jyc)
  ! enddo

  DO j=jdbeg,jdend
     DO i=idbeg,idend
        ux(i,j)=ux(i,j)+xicx(i)*(u(i,j)-u(i-1,j))/dxi
        uy(i,j)=uy(i,j)+etacy(j)*(uce(i,j)-uce(i,j-1))/deta

        vx(i,j)=vx(i,j)+xicx(i)*(vec(i+1,j)-vec(i-1,j))/twodxi
        vy(i,j)=vy(i,j)+etacy(j)*(v(i,j)-v(i,j-1))/deta

        ! d(i,j)=vx(i,j)
        d(i,j)=d(i,j)+xicx(i)*(u(i,j)-u(i-1,j))/dxi+etacy(j)*(v(i,j)-v(i,j-1))/deta
     ENDDO
  ENDDO

END SUBROUTINE uv_strain
