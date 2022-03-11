function X=proximalTV(x, lamTV);
% computes the solution to the minimisation problem
%
% min_{v}   1/2*||K(v)-x||_2^2+ lamTV*TV(v)
% with K(x)=x 
% via a primal dual hybrid method with extra gradient step for using the
% scheme for the indirect Tichonov problem from TVReg.m
%Prefereably use the direct method for computing prox_{TV}, which is
%implemented in proximalTV1.m or the implementation with DR-splitting
%proxTV_toolbox.m for better performance.
display=false;
sz=size(x);
x=reshape(x,sqrt(length(x(:))),[]);
N=size(x);
X=zeros(N); %initial primal variable
Xbar=X;

A= @(x) my_grad(x,1);         
Aadj= @(x,y) -my_div(x,y,1);  

Yx= zeros(N(1),N(1)); % initial TV component of dual variable
Yy= zeros(N(1),N(1));
R   =@(x) x;
Radj=@(x) x;
g=x;
Yz=zeros(size(g));   %define HS dual variable
Resolv_G=@(x) x ;
Resolv_F2=@(x,sigma) perform_soft_thresholding(x-sigma*g, sigma);                            
%set fix point iteration parameters
tau=1;
sigma1=1/(tau*8) ;
vec=normalize(randn(N));
sigma2=1/(tau*norm(R(vec))^2);

maxInter=1000;
epsilon=0.001;
cond=true; 
count=0;

while cond && count < maxInter
    [temp1, temp2]=A(Xbar);
    temp3=R(Xbar);
    
    %dual step
    Yx_temp=Yx+sigma1*temp1;
    Yy_temp=Yy+sigma1*temp2;
    Yz_temp=Yz+sigma2*temp3;
    [Yx_new,Yy_new]=Resolv_F1(Yx_temp,Yy_temp,lamTV); 
    Yz_new=Resolv_F2(Yz_temp,sigma2);                 
    % primal step
    X_temp=X-tau*Aadj(Yx_new,Yy_new)-tau*Radj(Yz_new);
    X_new=X_temp; 
    Xbar=2*X_new-X;  %extragradient
     
     if norm(X-X_new) < epsilon*norm(X) && norm(Yx-Yx_new)/norm(X) < epsilon ...      
          &&   norm(Yy-Yy_new)/norm(Yy) < eps && norm(Yz-Yz_new)/norm(Yz) < epsilon
         cond=false;
     end
    X=X_new;
    Yx=Yx_new;
    Yy=Yy_new;
    Yz=Yz_new;
    
     if mod(count,10)==0  && display  %display variables    
       montage({x/255,g/max(g,[],'all'),X/max(X,[],'all'),Yx/max(abs(Yx),[],'all'),Yy/max(Yy,[],'all'),Yz/max(Yz,[],'all')})
    end
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
            X=x;
            Y=y;
        else
        temp=max(nor,lamTV);
        X=lamTV.*x./temp;
        Y=lamTV.*y./temp;
    end
 
end