clear all
close all
clc

Vf = fopen('vobj100','w');

R=0.5;
n=32+1;
theta=2*pi/(n-1);
xs0=zeros(n,1);
ys0=zeros(n,1);
curv=zeros(n,1);

for i=1:n
    xs0(i)=cos((i-1)*theta)*R;
    ys0(i)=sin((i-1)*theta)*R;
    curv(i)=1/R;
end

plot(xs0,ys0);
axis([-1 1 -1 1]);
axis equal
fprintf(Vf,'&vertexinfo nvertex = 32, ndimension = 3, &end\n');
for i=1:n
    fprintf(Vf,'%16.12f %16.12f %16.12f \n',xs0(i),ys0(i),curv(i));
end

fclose(Vf);