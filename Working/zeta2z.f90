function zeta2z(zita) result(z)
double precision, intent(in):: zita
double precision:: z
double precision:: a,b,pi,e0,e1

pi=3.14159265358979323d0
e0=exp(0.0d0)
e1=exp(1.0d-2)

a=2.0d0*pi/(e1-e0)
b=pi*(e1+e0)/(e1-e0)

z = a*exp(zita)-b

end function zeta2z
