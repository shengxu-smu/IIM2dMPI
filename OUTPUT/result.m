clear all
close all

m=512;
n=m;
[p]=MPI_IO('p.dat',m,n);
[rhsp]=MPI_IO('rhsp.dat',m,n);

load xc.dat
load yc.dat

dx=2*pi/(m-1);
dy=2*pi/(m-1);

iop=pointinside(m,n);

p0=zeros(m,n);
for i=1:m
    for j=1:n
        if iop(i,j) ~= 0
            p0(i,j)=exp(-xc(i))*exp(-yc(j));
        else
            p0(i,j)=sin(xc(i))*sin(yc(j));
        end
    end
end
figure(1),
surf(p)
title('p')

error=abs(abs(p-p0));
figure(2),
surf(error)
title('error')
figure(3),
surf(p0)
title('p0')
figure(4),
surf(rhsp)
title('rhsp')
