MODULE stretchxy
  PRIVATE
  PUBLIC:: xi2x,x2xi,dxidx,dxdxi,d2xidx2,d2xdxi2,eta2y,y2eta,detady,dydeta,d2etady2,d2ydeta2
  DOUBLE PRECISION, PARAMETER:: pi=3.14159265358979323d0
  DOUBLE PRECISION, PARAMETER:: e =2.71828182845904524d0
  DOUBLE PRECISION, PARAMETER:: e0=1.0d0
  DOUBLE PRECISION, PARAMETER:: e1=3.0d0
  DOUBLE PRECISION:: a=2.0d0*pi/(e1-e0),b=pi*(e1+e0)/(e1-e0)

CONTAINS

  ! x(xi)
  FUNCTION xi2x(ksi) RESULT(x)
    DOUBLE PRECISION, INTENT(in):: ksi
    DOUBLE PRECISION:: x

    x=ksi
    !x = a*exp(ksi)-b
  END FUNCTION xi2x

  ! xi(x)
  FUNCTION x2xi(x) RESULT(ksi)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: ksi

    ksi=x
    !ksi=log((x+b)/a)
  END FUNCTION x2xi

  ! dxi/dx (xi)
  FUNCTION dxidx(ksi)
    DOUBLE PRECISION, INTENT(in):: ksi
    DOUBLE PRECISION:: dxidx

    dxidx=1
    !dxidx = exp(-ksi)/a
  END FUNCTION dxidx

  ! d2xi/dx2 (xi)
  FUNCTION d2xidx2(ksi)
    DOUBLE PRECISION, INTENT(in):: ksi
    DOUBLE PRECISION:: d2xidx2

    d2xidx2=0
    !d2xidx2 = -exp(2.0d0*ksi)/(a*a)
  END FUNCTION d2xidx2

  ! dx/dxi (x)
  FUNCTION dxdxi(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: dxdxi

    dxdxi=1
    !dxidx = exp(-ksi)/a
  END FUNCTION dxdxi

  ! d2x/dxi2 (x)
  FUNCTION d2xdxi2(x)
    DOUBLE PRECISION, INTENT(in):: x
    DOUBLE PRECISION:: d2xdxi2

    d2xdxi2=0
    !d2xidx2 = -exp(2.0d0*ksi)/(a*a)
  END FUNCTION d2xdxi2

  ! y(eta)
  FUNCTION eta2y(ita) RESULT(y)
    DOUBLE PRECISION, INTENT(in):: ita
    DOUBLE PRECISION:: y

    y=ita
    !y = a*exp(ita)-b
  END FUNCTION eta2y

  ! eta(y)
  FUNCTION y2eta(y) RESULT(ita)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: ita

    ita=y
    !ita=log((y+b)/a)
  END FUNCTION y2eta

  ! deta/dy (eta)
  FUNCTION detady(ita)
    DOUBLE PRECISION, INTENT(in):: ita
    DOUBLE PRECISION:: detady

    detady=1
    !detady = exp(-ita)/a
  END FUNCTION detady

  ! d2eta/dy2 (eta)
  FUNCTION d2etady2(ita)
    DOUBLE PRECISION, INTENT(in):: ita
    DOUBLE PRECISION:: d2etady2

    d2etady2=0
    !d2etady2 = -exp(2.0d0*ita)/(a*a)
  END FUNCTION d2etady2

  ! dy/deta (y)
  FUNCTION dydeta(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: dydeta

    dydeta=1
    !detady = exp(-ita)/a
  END FUNCTION dydeta

  ! d2y/deta2 (y)
  FUNCTION d2ydeta2(y)
    DOUBLE PRECISION, INTENT(in):: y
    DOUBLE PRECISION:: d2ydeta2

    d2ydeta2=0
    !d2etady2 = -exp(2.0d0*ita)/(a*a)
  END FUNCTION d2ydeta2

END MODULE stretchxy
