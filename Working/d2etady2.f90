function d2etady2(ita)
double precision, intent(in):: ita
double precision:: d2etady2
double precision:: a,b,pi,e0,e1

pi=3.14159265358979323d0
e0=exp(0.0d0)
e1=exp(1.0d-2)

a=2.0d0*pi/(e1-e0)
b=pi*(e1+e0)/(e1-e0)

d2etady2 = -exp(2.0d0*ita)/(a*a)

end function d2etady2
