function eta2y(ita) result(y)
double precision, intent(in):: ita
double precision:: y
double precision:: a,b,pi,e0,e1

pi=3.14159265358979323d0
e0=exp(0.0d0)
e1=exp(1.0d-2)

a=2.0d0*pi/(e1-e0)
b=pi*(e1+e0)/(e1-e0)

y = a*exp(ita)-b

end function eta2y
