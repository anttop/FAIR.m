function Div = div(Fx,Fy,h)
%compute divergence via a backwards diffrences
%with zeros values at the boundary (constant at the boundary)

%h stepsize
%m=size(F);
%Fx=reshape(F(1,:,:),m(2),m(3));
%Fy=reshape(F(2,:,:),m(2),m(3));
help1=circshift(Fx,1,1);
Div1=(Fx-help1)/h;
Div1(1,:)=Fx(1,:)/h;
Div1(end,:)=-help1(end,:)/h;

help2=circshift(Fy,1,2);
Div2=(Fy-help2)/h;
Div2(:,1)=Fy(:,1);
Div2(:,end)=-help2(:,end)/h;

Div=Div1+Div2;

