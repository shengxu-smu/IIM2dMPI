subroutine Cartesiangrid
use para
use Euler
integer:: i,j,k
double precision:: xi0,xi1,ksi,eta0,eta1,ita,zeta0,zeta1,zita

interface
  function x2xi(x) result(ksi)
    double precision, intent(in):: x
    double precision:: ksi
  end function x2xi
  function y2eta(y) result(ita)
    double precision, intent(in):: y
    double precision:: ita
  end function y2eta
  function z2zeta(z) result(zita)
    double precision, intent(in):: z
    double precision:: zita
  end function z2zeta

  function xi2x(ksi) result(x)
    double precision, intent(in):: ksi
    double precision:: x
  end function xi2x
  function eta2y(ita) result(y)
    double precision, intent(in):: ita
    double precision:: y
  end function eta2y
  function zeta2z(zita) result(z)
    double precision, intent(in):: zita
    double precision:: z
  end function zeta2z

  function dxidx(ksi)
    double precision, intent(in):: ksi
    double precision:: dxidx
  end function dxidx
  function detady(ita)
    double precision, intent(in):: ita
    double precision:: detady
  end function detady
  function dzetadz(zita)
    double precision, intent(in):: zita
    double precision:: dzetadz
  end function dzetadz

  function d2xidx2(ksi)
    double precision, intent(in):: ksi
    double precision:: d2xidx2
  end function d2xidx2
  function d2etady2(ita)
    double precision, intent(in):: ita
    double precision:: d2etady2
  end function d2etady2
  function d2zetadz2(zita)
    double precision, intent(in):: zita
    double precision:: d2zetadz2
  end function d2zetadz2
end interface

! x=x(xi)
xi0=x2xi(x0)
xi1=x2xi(x1)
!dxi=(xi1-xi0)/dble(nx-1)
dxi=(xi1-xi0)/dble(nx+1) ! for testing hypre w. Dirichlet BCs
halfdxi=HALF*dxi
twodxi=TWO*dxi
do i=ipbeg,ipend
! ksi=xi0+(dble(i)-ONE)*dxi  
  ksi=xi0+dble(i)*dxi ! for testing hypre w. Dirichlet BCs
  xic(i)=ksi
  xicx(i)=dxidx(ksi)
  xicxx(i)=d2xidx2(ksi)
  xc(i)=xi2x(ksi)
enddo
do i=iubeg,iuend
  ksi=xi0+(dble(i)-HALF)*dxi
  xif(i)=ksi
  xifx(i)=dxidx(ksi)
  xifxx(i)=d2xidx2(ksi)
  xf(i)=xi2x(ksi)
enddo

! y=y(eta)
eta0=y2eta(y0)
eta1=y2eta(y1)
!deta=(eta1-eta0)/dble(ny-1)
deta=(eta1-eta0)/dble(ny+1) ! for testing hypre w. Dirichlet BCs
halfdeta=HALF*deta
twodeta=TWO*deta
do j=jpbeg,jpend
! ita=eta0+(dble(j)-ONE)*deta
  ita=eta0+dble(j)*deta ! for testing hypre w. Dirichlet BCs
  etac(j)=ita
  etacy(j)=detady(ita)
  etacyy(j)=d2etady2(ita)
  yc(j)=eta2y(ita)
enddo
do j=jvbeg,jvend
  ita=eta0+(dble(j)-HALF)*deta
  etaf(j)=ita
  etafy(j)=detady(ita)
  etafyy(j)=d2etady2(ita)
  yf(j)=eta2y(ita)
enddo

! z=z(zeta)
zeta0=z2zeta(z0)
zeta1=z2zeta(z1)
!dzeta=(zeta1-zeta0)/dble(nz-1)
dzeta=(zeta1-zeta0)/dble(nz+1)! for testing hypre w. Dirichlet BCs
halfdzeta=HALF*dzeta
twodzeta=TWO*dzeta
do k=kpbeg,kpend
! zita=zeta0+(dble(k)-ONE)*dzeta
  zita=zeta0+dble(k)*dzeta! for testing hypre w. Dirichlet BCs
  zetac=zita
  zetacz(k)=dzetadz(zita)
  zetaczz(k)=d2zetadz2(zita)
  zc(k)=zeta2z(zita)
enddo
do k=kwbeg,kwend
  zita=zeta0+(dble(k)-HALF)*dzeta
  zetaf=zita
  zetafz(k)=dzetadz(zita)
  zetafzz(k)=d2zetadz2(zita)
  zf(k)=zeta2z(zita)
enddo

r2dxideta=(dxi*dxi)/(deta*deta)
r2dxidzeta=(dxi*dxi)/(dzeta*dzeta)

end subroutine Cartesiangrid
