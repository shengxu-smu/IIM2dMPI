clear all
close all
clc

Vf = fopen('vobj25','w');

a=1;
ns=1024;
l=4*a;
n=ns+1;dn=(n-1)/8;

xs0=zeros(n,1);
ys0=zeros(n,1);
curv=zeros(n,1);
alfa=zeros(n,1);
dx=l/ns;

for i=1:dn+1
    xs0(i)=a/2;
    ys0(i)=(i-1)*dx;
end
for i=dn+2:3*dn+1
    xs0(i)=a/2-(i-dn-1)*dx;
    ys0(i)=a/2;
end
for i=3*dn+2:5*dn+1
    xs0(i)=-a/2;
    ys0(i)=a/2-(i-3*dn-1)*dx;
end
for i=5*dn+2:7*dn+1
    xs0(i)=-a/2+(i-5*dn-1)*dx;
    ys0(i)=-a/2;
end
for i=7*dn+2:n
    xs0(i)=a/2;
    ys0(i)=-a/2+(i-7*dn-1)*dx;
end

plot(xs0,ys0);
axis([-1 1 -1 1])
axis equal
fprintf(Vf,'&vertexinfo nvertex = 1024, ndimension = 3, &end\n');
for i=1:n
    fprintf(Vf,'%16.12f %16.12f %16.12f\n',xs0(i),ys0(i),curv(i));
end

fclose(Vf);
