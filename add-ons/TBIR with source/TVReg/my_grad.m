function [Fx, Fy] = grad(F,h);
%computes discrete gradient of F with a forward diffrence, 
%where the value is set to zero on the right boundary
%h determines the stepsize

Fx=(circshift(F,-1,1)-F)/h;
Fx(:,end)=0;

Fy=(circshift(F,-1,2)-F)/h;
Fy(end,:)=0;

m=size(F);
Fout=zeros(2,m(1), m(2));
FoutFout(1,:,:)=Fx;
Fout(2,:,:)=Fy;

