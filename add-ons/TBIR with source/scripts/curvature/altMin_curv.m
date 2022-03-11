%This script solves the indirect image registration problem  with source
%term (TV regularisation) via an
%alternating minimisation scheme over the variables (v,z),
%where a curvature regulariser is posed on v.
%For every iteration:
% i)  An inexact Gauss-Newton step is performed on v for fixed z.
% ii) The variable v is fixed and the indirect Tikhinov problem with TV
% regularisation over z is solved via a primal dual hybrid-gradient ( see
% TVReg.m)

clear;
close all;
clc;
plot=false;
loadData = false;
name='test';
% Load images.
path = fullfile(FAIRpath, 'add-ons', 'TBIR with source', 'data');
file1 = 'phantom.png';
file2 = 'deform3.png';
file3 = 'deformphantom.png';
image1 = double(imresize(imread(fullfile(path, file1)), [128, 128]));
image2 = double(imresize(imread(fullfile(path, file3)), [128, 128]));
image3 = double(imresize(imread(fullfile(path, file3)), [128, 128]));
m = size(image1);
% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'TBIR with source', 'results', 'curvature',  'alternating minimisation');
mkdir(outputfolder);
% Save template, unknown image, and measurements to results folder
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(image2 / 255, fullfile(outputfolder, sprintf('%s_target.png', name)));


% Name of dataset.
name ='phantom';
% Save size of template.
m = size(image1);
% Set bumber of Runge-Kutta steps.
N = 5;
% Define data term.
dist = 'SSD_op';
% Define regularization term.
reg = 'mfCurvatureST';
% Define objective.
objfun = 'LDDMMobjFctn';
% Define image model.
imageModel = 'splineInterMex';
% Define temporal discretization of the velocity (number of time steps is
% then nt + 1.
nt = 1;
% Define noise level.
noise = 0.0;
% Set regularization parameters.
alpha = [1,5];

lamTV=5;
% Set Hessian shift.
hessianShift = 1e-2;
% Set domain size.
omega = [0, 1, 0, 1];
% Set domain for velocities by padding.
pad = 0.5;
omegaV = omega;
omegaV(1:2:end) = omegaV(1:2:end) - pad;
omegaV(2:2:end) = omega(2:2:end) + pad;

% Initialize models.
imgModel('reset', 'imgModel', imageModel)
trafo('reset', 'trafo', 'affine2D')
distance('reset', 'distance', dist);
viewImage('reset', 'viewImage', 'viewImage2D', 'colormap', gray(256));
NPIRpara = optPara('NPIR-GN');
NPIRpara.maxIter = 2;
NPIRpara.scheme = @GaussNewtonLDDMM;

% Create multilevel versions of template.
[ML, ~, maxLevel, ~] = getMultilevel({image1}, omega, m, 'fig', 0);
minLevel=4;
% Set directions for Radon transform.
theta = linspace(10, 180, 10);

% Set up operators for all levels.
for k=minLevel:maxLevel
    [ML{k}.K, ML{k}.Kadj, ML{k}.cleanup, ML{k}.ndet] = createRadon2d(size(ML{k}.T), theta, 2^(-k + minLevel));
end

% Run algorithm for each setting.
mV = @(m) ceil(1*m);

% Check if measurements exists, otherwise create.
sinogramfile = fullfile(outputfolder, 'Sinograms', sprintf('%s_measurements_%g.mat', name, noise));
if(exist(sinogramfile, 'file'))
    S = load(sinogramfile);
    ML{maxLevel}.R = S.R;
else
    % Apply operator on finest level to generate synthetic measurements.
    R = ML{maxLevel}.K(image2);
    
    % Add noise to measurements.
    R = addnoise(R, noise);
    ML{maxLevel}.R = R;

    % Save measurements.
    mkdir(fullfile(outputfolder, 'Sinograms'));
    save(sinogramfile, 'R', 'theta', 'm');
end

% Check if measurements exists, otherwise create.
sinogramfile = fullfile(outputfolder, 'Sinograms', sprintf('%s_sino_%g_%g.mat', name, noise, length(theta)));
if(exist(sinogramfile, 'file'))
    S = load(sinogramfile);
    ML{maxLevel}.R = S.R;
else
    % Apply operator on finest level to generate synthetic measurements.
    R = ML{maxLevel}.K(image2);

    % Add noise to measurements.
    R = addnoise(R, noise);
    ML{maxLevel}.R = R;

    % Save measurements.
    mkdir(fullfile(outputfolder, 'Sinograms'));
    save(sinogramfile, 'R', 'theta', 'm');
end

%save measurements
Rsize = size(ML{maxLevel}.R, 2);
Rsq = imresize(ML{maxLevel}.R, [Rsize, Rsize], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino_%.2f_%f.png', name, noise,length(theta))));


% Create multilevel versions of measurements.
ML = multilevelRadon2d(ML, maxLevel, minLevel);

vcc=[];

ticId=tic;
for k=minLevel:maxLevel
    %set TV-Reg parameters
    X=zeros(ML{k}.m); %initial primal variable
    Xbar=X;
    R=ML{k}.R;
    
    A= @(x)my_grad(x,1);         
    Aadj= @(x,y) -my_div(x,y,1); 
    Yx= zeros(ML{k}.m); % initial TV component of dual variable
    Yy= zeros(ML{k}.m);
    
    %set forward operator for TV-step
    R1=ML{k}.K;
    R1adj=@(x) reshape(ML{k}.Kadj(x),ML{k}.m);

    g=ML{k}.R-R1(ML{k}.T);
    R=ML{k}.R;
    Yz=zeros(size(g));
 
    Resolv_G=@(x) x ;
    Resolv_F2=@(x,sigma) perform_soft_thresholding(x, sigma)

    %set fix point iteration parameters
    tau=1;
    sigma1=1/(tau*8) ;
    vec=normalize(randn(ML{k}.m));
    sigma2=1/(tau*norm(R1(vec))^2);

    epsilon=0.001; %error tolerance
    STOP=false;
    i=0;
    condTV=true; 
    maxInterTV=100;
   
    %interpolate velocity form last level
    if k~=minLevel
        vnew=reshape(vc,[],2);
        v0=[mfPu(vnew(:,1),omega,mV(ML{k}.m/2)) mfPu(vnew(:,2),omega,mV(ML{k}.m/2))];
        vcc=v0(:);
    else
       vcc= getVelocityStartingGuess(omegaV,mV(ML{k}.m),nt);
    end
    rec=[];
    maxInter=1;
    
    while  (i<maxInter) && ~STOP
        v_old=vcc;
        X_old=X;
        countTV=0;
        
        %run TV step
        while condTV & countTV < maxInterTV

        [temp1, temp2]=A(Xbar);
        temp3=R1(Xbar)-g;
        %dual step
        Yx_temp=Yx+sigma1*temp1;
        Yy_temp=Yy+sigma1*temp2;
        Yz_temp=Yz+sigma2*temp3;   
        [Yx_new,Yy_new]=Resolv_F1(Yx_temp,Yy_temp,lamTV);
        Yz_new=Resolv_F2(Yz_temp,sigma2);


        X_temp=X-tau*Aadj(Yx_new,Yy_new)-tau*R1adj(Yz_new);
        X_new=X_temp; % primal step
        Xbar=2*X_new-X;  %extragradient

        if norm(X-X_new)/norm(X) < epsilon & norm(Yx-Yx_new)/norm(X) < epsilon ...      
                norm(Yy-Yy_new)/norm(Yy) < eps & norm(Yz-Yz_new)/norm(Yz) < epsilon
            cond=false;
            break
        end
        X=X_new;
        Yx=Yx_new;
        Yy=Yy_new;
        Yz=Yz_new;

          if mod(countTV,10)==0  && plot==true; %display 
            montage({rec/255,g/max(g,[],'all'),X/max(X,[],'all'),Yx/max(abs(Yx),[],'all'),Yy/max(Yy,[],'all'),Yz/max(Yz,[],'all')})
         end
       countTV=countTV+1;
       
     end         
    ML{k}.R=R-R1(X);        
%      Run indirect registration.
    regularizer('reset', 'regularizer', reg, 'nt', nt,...
        'alpha', alpha, 'HessianShift', hessianShift);
    [vc, ~, ~, his] = MLLDDMM(ML, 'operator', true, 'minLevel',...
        k, 'maxLevel', k, 'omegaV', omegaV, 'mV', mV,...
        'N', N, 'parametric', false, 'NPIRpara', NPIRpara,...
        'NPIRobj', str2func(objfun), 'plots', plot,'vRef', vcc);
    
     STOP = (norm(vcc-vc) < epsilon*norm(vc)) && (norm(X-X_old) < epsilon*norm(X));
     vcc=vc;
%      Transform template and reshape.
    yc = getTrafoFromInstationaryVelocityRK4(vc, getNodalGrid(omega,ML{k}.m),...
        'omega', omegaV, 'm', ML{k}.m, 'nt', nt, 'tspan', [1, 0], 'N', N);
    rec = linearInterMex(ML{maxLevel}.T, omega, center(yc, ML{k}.m));
    rec = reshape(rec, ML{k}.m);

    g=R-R1(rec);
    
    i=i+1;
    end
   
end
      
elapsed=toc(ticId);
%  Transform template and reshape.
yc = getTrafoFromInstationaryVelocityRK4(vc, getNodalGrid(omega,m),...
    'omega', omegaV, 'm', m, 'nt', nt, 'tspan', [1, 0], 'N', N);
rec = linearInterMex(ML{maxLevel}.T, omega, center(yc, m));
rec = reshape(rec, m);
z=reshape(X,m);
% Output stats.
fprintf('Elapsed time is: %.2f seconds, SSIM=%.3f.\n',elapsed, ssim(rec, image2));

%save results
[resfile,resfile1,resfile2,resfile3,resfile4,paramfile] = saveres(name, outputfolder, image1,...
    image2,image3, rec,z, dist, reg, objfun, imageModel, N, nt, alpha, theta,...
    noise, elapsed, 0, lamTV);

% Free resources.
for k=minLevel:maxLevel
    ML{k}.cleanup();
end




function [X,Y]=Resolv_F1(x,y,lamTV)
        nor=(x.^2.+y.^2);
        nor=sqrt(nor);
        if nor==0 ;
            X=x;
            Y=y;
        elseif lamTV==0;
            X=0.*x;
            Y=0.*y;
        else
        temp=max(nor/lamTV,1);
        X=x./temp;
        Y=y./temp;
    end
    
end