function [iop]=pointinside(m,n)
close all

fid=fopen('iop.dat');
iop=fread(fid,[m n],'int16');

% fid=fopen('iou.dat');
% iou=fread(fid,[m n],'int16');
% 
% fid=fopen('iov.dat');
% iov=fread(fid,[m n],'int16');
% fclose('all');
% 
load object1.dat
load object2.dat
load xc.dat
load yc.dat
load xf.dat
load yf.dat

object=object1;
% iop
% figure,
% hold on
% plot(object(:,1),object(:,2),'k-')
% for i=1:m
%     for j=1:n
%         if(iop(i,j)==1)
%             plot(xc(i),yc(j),'.r');
%         end
%     end
% end
% hold off
% title('iop');
% xlabel('x');ylabel('y');
% % axis([-1 1 -1 1])
% axis equal

% iou
% figure,
% hold on
% plot(object2(:,1),object2(:,2),'k-')
% for i=1:m
%     for j=1:n
%         if(iou(i,j)==1)
%             plot(xf(i+1),yc(j),'.r');
%         end
%     end
% end
% hold off
% % axis([-1 1 -1 1])
% title('iou');
% xlabel('x');ylabel('y');
% axis equal
% 
% % iov
% figure,
% hold on
% plot(object2(:,1),object2(:,2),'k-')
% for i=1:m
%     for j=1:n
%         if(iov(i,j)==1)
%             plot(xc(i),yf(j+1),'.r');
%         end
%     end
% end
% hold off
% % axis([-0.6 0.6 -0.6 0.6])
% title('iov');
% xlabel('x');ylabel('y');
% axis equal

end