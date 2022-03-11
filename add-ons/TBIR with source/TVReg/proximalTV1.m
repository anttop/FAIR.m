function X=proximalTV1(x, lamTV);
% This script solves the problem 
%
% min_{} 1/2 ||x-g||_2^2 +lamTV*TV(x)
% via a primal dual hybrid gradient method with extragradient.

sz=size(x); %save size of input variable
x=reshape(x,sqrt(length(x(:))),[]);
N=size(x);
X=zeros(N); %initial primal variable
Xbar=X;

A= @(x) my_grad(x,1);         
Aadj= @(x,y) -my_div(x,y,1);  

Yx= zeros(N(1),N(1)); % initial dual variable
Yy= zeros(N(1),N(1));
Yz=zeros(size(x));   %define HS dual variable
%set fix point iteration parameters
sigma1=1;
tau=1/8;

maxInter=1000;
epsilon=0.001;
cond=true; 
count=0;

while cond & count < maxInter
    [temp1, temp2]=A(Xbar);
    %dual step
    Yx_temp=Yx+sigma1*temp1;
    Yy_temp=Yy+sigma1*temp2;
    [Yx_new,Yy_new]=Resolv_F1(Yx_temp,Yy_temp,lamTV);
    % primal step
    X_temp=X-tau*Aadj(Yx_new,Yy_new);
    X_new= (X_temp+tau*x)./(1+tau); 
    Xbar=2*X_new-X;  %extragradient
     
    if norm(X-X_new) < epsilon*norm(X) && norm(Yx-Yx_new) < epsilon*norm(Yx) && norm(Yy-Yy_new) < eps*norm(Yy)
        cond=false;
    end
    X=X_new;
    Yx=Yx_new;
    Yy=Yy_new;

    
%       if mod(count,10)==0  %display variables
%        montage({g/max(g,[],'all'),X/max(X,[],'all'),Yx/max(abs(Yx),[],'all'),Yy/max(Yy,[],'all')})
%      end
    count=count+1;
end
 X=reshape(X,sz); 
end

function [X,Y]=Resolv_F1(x,y,lamTV)
        if lamTV==0;
            X=0.*x;
            Y=0.*y;
            return
        end
        nor=(x.^2.+y.^2);
        nor=sqrt(nor);
        if nor==0 ;
            X=lamTV.*x;
            Y=lamTV.*y;
        else
        temp=max(nor,lamTV);
        X=lamTV*x./temp;
        Y=lamTV*y./temp;
    end
    
end

function a=soft_thresholding(x,sigma)
    nor=x.^2
    nor=sqrt(nor);
    if nor==0;
        a=0.*x;
    else
     temp=max(0,nor-sigma);
     temp=temp./nor;
     a=temp.*x;
    end
end