SUBROUTINE FiniteDiff(xx,yy,yfirst,ysecond)

  USE mpi
  USE para
  USE Euler

  DOUBLE PRECISION, DIMENSION(0:3)::xx,yy
  DOUBLE PRECISION:: yfirst,ysecond
  DOUBLE PRECISION,DIMENSION(1:2)::a,b,c,cof1,cof2,cof3,cof0,dd

  a=ZERO
  b=ZERO
  c=ZERO
  dd=ZERO
  cof1=ZERO
  cof2=ZERO
  cof3=ZERO
  cof0=ZERO

  a(1)=-(xx(2)-xx(0))*(xx(2)-xx(0))/(xx(1)-xx(0))/(xx(1)-xx(0))
  b(1)=1.0d0
  dd(1)=(a(1)*(xx(1)-xx(0))+b(1)*(xx(2)-xx(0)))
  cof1(1)=a(1)/dd(1)
  cof2(1)=b(1)/dd(1)
  cof0(1)=-(a(1)+b(1))/dd(1)

  a(2)=(xx(2)-xx(0))*(xx(3)+xx(2)-2*xx(0))*(xx(3)-xx(2))/ &
       (xx(1)-xx(0))/(xx(3)+xx(1)-2*xx(0))/(xx(1)-xx(3))
  b(2)=1.0d0
  c(2)=-a(2)*(xx(0+1)-xx(0))/(xx(0+3)-xx(0)) &
       -b(2)*(xx(0+2)-xx(0))/(xx(0+3)-xx(0))
  dd(2)=0.5d0*(a(2)*(xx(1)-xx(0))*(xx(1)-xx(0))+&
       b(2)*(xx(2)-xx(0))*(xx(2)-xx(0))+&
       c(2)*(xx(3)-xx(0))*(xx(3)-xx(0)))
  cof1(2)=a(2)/dd(2)
  cof2(2)=b(2)/dd(2)
  cof3(2)=c(2)/dd(2)
  cof0(2)=-(a(2)+b(2)+c(2))/dd(2)

  yfirst  = cof0(1)*yy(0)+cof1(1)*yy(1)+cof2(1)*yy(2)
  ysecond = cof0(2)*yy(0)+cof1(2)*yy(1)+cof2(2)*yy(2)+cof3(2)*yy(3)

END SUBROUTINE FiniteDiff
