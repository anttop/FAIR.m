% This script solves the proposed problem with a curvature regulariser instead of
% third order regularisation.
%__________________________________________________________________________
noise=0; %noise level
loadData = true;
plot=true;
name='lotus';
% Load images.
path = fullfile(FAIRpath, 'add-ons', 'TBIR with source','data');
% Load measurements.
D = load('sinogram.mat');
sinogram = D.sinogram;
% Correct for centre of rotation and downsize.
sinogram = sinogram(1:2221,1:360);
factor = 128 / 1500;
sinogram = imresize(sinogram, [ceil(factor * 2221), 360], 'bilinear', 'Antialiasing', false);

% Reconstruct image using FBP as ground truth.
D = 2240 * 54 / 12;
FSS = 54 / 63;
G = 2240 * (63 - 54)/12;
image2 = ifanbeam(sinogram, D,...       
    'FanRotationIncrement', 1,...
    'FanSensorGeometry', 'line',...
    'FanSensorSpacing', FSS,...
    'OutputSize', 128);              %target image
m = size(image2);

% The template image is based on the following deformation of the target
% image image2. The additional structures were added manuanlly with GIMP.
% The XCF-file can be found in FAIR/add-ons/experiments/data
% % Introduce a deformation to the template.
% sigma = 20;
% c = [70, 45];
% [X, Y] = ndgrid(1:m(1), 1:m(2));
% K = exp(-((X - c(1)).^2 + (Y - c(2)).^2) / (2 * sigma^2));
% v = 3 * cat(3, -4 * K, 1 * K);
% image3 = imwarp(image2, v);   %deformation of the template

%load template
file1 = 'lotus_source2.png';
image1 = double(imresize(imread(fullfile(path, file1)), [128, 128]))./255;
image1= max(image2(:))*image1;
m = size(image1);

%reduce data form sinogram
sinogram = sinogram(:, 1:15:180);
theta = 1:15:180;
% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'TBIR with source', 'results','curvature','lotus');
mkdir(outputfolder);
% Save template, unknown image, and measurements to results folder
imwrite(uint8(255 * image1 / max(image1(:))), fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(uint8(255 * image2 / max(image2(:))), fullfile(outputfolder, sprintf('%s_target.png', name)));

% Set parameters and models
N = 5; % Runge-Kutta steps
nt = 1; % Temporal discretization of velocity
noise = 0.0; % Noise
alpha = [1,1]; % Regularization velocity.
lambda = 0.001; % Regularization source
dist = 'NCC_op';
reg = 'mfCurvatureST';
objfun = 'PALM_NCC_objFctn';
imageModel = 'splineInterMex';

% Create domains
omega = [0, 1, 0, 1];
pad = 0.5;
omegaV = omega;
omegaV(1:2:end) = omegaV(1:2:end) - pad;
omegaV(2:2:end) = omega(2:2:end) + pad;
mV = @(m) ceil(1*m);

% Initialize models.
imgModel('reset', 'imgModel', imageModel)
trafo('reset', 'trafo', 'affine2D')
distance('reset', 'distance', dist);
viewImage('reset', 'viewImage', 'viewImage2D', 'colormap', gray(256));
regularizer('reset', 'regularizer', reg, 'nt', nt,...
    'alpha', alpha, 'HessianShift', 0);
shift = 1e-3;

% Create multilevel versions of template.
[ML, ~, maxLevel, ~] = getMultilevel(image1, omega, m, 'fig', 0);
% Set starting level.
minLevel = 6;
% Set up operators for all levels.
ndet = size(sinogram, 1);
for k=minLevel:maxLevel
    [ML{k}.K, ML{k}.Kadj, ML{k}.cleanup, ML{k}.ndet] = createRadon2d_fan(size(ML{k}.T), theta, ceil(ndet * 2^(k - maxLevel)), 2^(-k + minLevel), D, G);
end

% Run algorithm for each setting.
mV = @(m) ceil(1*m);

% Set measurements.
ML{maxLevel}.R = permute(sinogram, [2, 1]);

% Save measurements to results folder.
Rsize = size(ML{maxLevel}.R, 2);
Rsq = imresize(ML{maxLevel}.R, [Rsize, Rsize], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino_%.2f.png', name, noise)));

% Create multilevel versions of measurements.
ML = multilevelRadon2d_imresize(ML, maxLevel, minLevel);




%PALM
ProxG = @(x,sigma) proximalTV1(x,sigma*lambda); %proximal map for second variable
maxIter=100;
tol=10^-3;
vRef=[];
eta=3; % Backtracking parameter

ticId = tic; 
for level=minLevel:maxLevel
    
    % Update m,grid,data and coefficients
    m    = ML{level}.m;
    xc   = getCellCenteredGrid(omega,m);
   [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
    Rc = ML{level}.R;
    K=ML{level}.K;
    Kadj=ML{level}.Kadj;
    trafo([],xc)
    
    % compute starting guess y0 for this level and stopping value yStop
    if level == minLevel
        if isempty(vRef)
            vRef = getVelocityStartingGuess(omegaV,mV(m),nt);
        end
        v0 = vRef;  % the best known so far
        z = reshape(zeros(m),[],1);
    else
        % prolongate vc (coarse) vRef (current)
        vcOld   = reshape(vc,[],nt+1);
        vRefOld = reshape(vRef,[],nt+1);
        v0 = []; vRef = [];
        for k=1:nt+1
            v0   = [v0 mfPu(vcOld(:,k),omega,mV(m/2))];
            vRef = [vRef mfPu(vRefOld(:,k),omega,mV(m/2))];
        end
        v0 = v0(:); vRef = vRef(:);
        z = imresize(reshape(z,m./2),m,'nearest');
        z=reshape(z,[],1);
    end
    vStop = 0*getVelocityStartingGuess(omegaV,mV(m),nt);
    
    % Initialize iPALM run
    v=v0; % Starting guess
    vi=v;
    zi=z;
    J=0;
    Jz=0;
    iter=1;
    Knorm=operator_norm(K,Kadj,rand(m(1)*m(2),1));
    sigma=Knorm^2; %step size second variable
    tau=sigma;     %step size first variable
    sigma=1;
    %proximal mapping for first variable
    [~,~,ScdDev] = regularizer(v0,omegaV,m,'doDerivative',1,'tspan',[1 0]);
    ProxF = @(x)x + ScdDev.d2S(x,omega,m)/tau;
    Preconditioner = @(x) solveSpectral(omega,m,[1 0],nt,1,1/tau,x);
    STOP(1) = false;
    STOP(2) = false;
    count=0;
    count1=0;
    Amplitude = @(u) sqrt(sum(u.^2, 3));
    while ~any(STOP)
        % Keep record
        v_old = v;
        z_old = z;
        J_old = J;
        Jz_old=Jz;
        
        % Update v with line search
        NPIRfctn = @(vc) PALM_NCC_objFctn(T,Rc,z,omega,m,vRef,xc,omegaV,mV(m),N,vc,K,Kadj);
        J=NPIRfctn(v);
         %determine stepsize tau via backtracking scheme
        [~,~,dJ,~]=NPIRfctn(vi);
        Jhat = max(J,J_old);
        BT_cond=false;
        while not(BT_cond)
            v_temp = vi-dJ'./tau;
            v_temp = pcg(ProxF,v_temp,10^-5,300,Preconditioner,[],vi);
            w=v_temp-v;
            J1=NPIRfctn(v_temp);
            if J1 <= Jhat+dJ*w+(tau/2)*(norm(w))^2
                BT_cond= true;
                count=count+1;
                % Sometimes try to increase tau
                if count == 5
                    count = 0;
                    tau = tau/eta;
                    ProxF = @(x)x + ScdDev.d2S(x,omega,m)/tau;
                    Preconditioner = @(x) solveSpectral(omega,m,[1 0],nt,1,1/tau,x);
                end
            else
                tau=eta*tau;
                ProxF = @(x)x + ScdDev.d2S(x,omega,m)/tau;
                Preconditioner = @(x) solveSpectral(omega,m,[1 0],nt,1,1/tau,x);
            end
        end
        J1+lambda * sum(sum(Amplitude(grad(z)))); %monitor function value
        v = v_temp;
        vi = v + (iter-1)/(iter+2).*(v-v_old);   %inertia
        
        % compute deformation part
        yc = getTrafoFromInstationaryVelocityRK4(v, getNodalGrid(omega,m),...
            'omega', omegaV, 'm', m, 'nt', nt, 'tspan', [1 0], 'N', N);
        [rec,~] = imgModel(T, omega, center(yc, m));
        
       
        % second variable (backtracking via line search
         eta1=1.3;
         NPIRfctn = @(z) PALM_NCC_objFctn(T,Rc,z,omega,m,vRef,xc,omegaV,mV(m),N,v,K,Kadj);
         Jz=NPIRfctn(z)/1e6; %function value 
         f=K(rec+zi);
         g=Rc/norm(Rc, 'fro');
         f_sqnr=(f(:)'*f(:));
         scprod= f(:)'*g(:);
         dD=2*Kadj((-scprod/f_sqnr)*g+f*(scprod/f_sqnr)^2); %gradient
          
         Jhat_z = max(Jz,Jz_old);
         BT_cond1=false;
          while not(BT_cond1)
                
            z_temp = zi-dD/sigma;
            z_temp = ProxG(z_temp,sigma);
            x=z_temp-z;
            Jz1=NPIRfctn(z_temp)/1e6;
                if Jz1 <= Jhat_z+dD'*x+(sigma/2)*(norm(x))^2
                    BT_cond1= true;
                    count1=count1+1;
                    % Sometimes try to increase sigma
                    if count1 == 5
                        count1 = 0;
                        sigma = sigma/eta1;  
                    end
                else
                    sigma=eta1*sigma;
                    if (sigma>1e20)
                        z_temp=z;
                        BT_cond1=true;
                    end        
                end
            end
          z=z_temp;
          zi = z + (iter-1)/(iter+2).*(z-z_old); %inertia
       

        % Create plot
        A=reshape(rec, m);           %deformation part
        AA=reshape(rec+z, m);        %reconstruction
       % 'maxdef'; max(rec); 'max z';max(z(:)); 'max rec'; max(AA(:));
        if plot && mod(iter,1)==0 
          %display images
          montage({rescale(A),reshape(z,m)/(max(abs(z(:)))), reshape(-z,m)/(max(abs(z(:)))), rescale(AA)});
        end
        % Set variables for next iteration
        iter = iter + 1;
        STOP(1) = iter > maxIter;
        STOP(2) = (norm(v-v_old) < tol*norm(v)) && (norm(z-z_old) <= tol*norm(z));
     %   STOP(3) = (abs(J-J_old) < tol*abs(J));
    end
    
    %postprocessing: inexact Gauss-Newton on variable v 
    % Define objective.
    objective_new = @(vc) NCC_objFctn(T,Rc,z,omega,m,vRef,xc,omegaV,mV(m),N,vc,K,Kadj);
    for ii=1:15
        [J,~,dJ,H] = objective_new(v);
        J + lambda * sum(sum(Amplitude(grad(z))))
        P       = H.d2D.P; % grid-to-grid operator
        dr      = H.d2D.dr;
        dr_adj  = H.d2D.dr_adj;
        d2psi   = H.d2D.d2psi;
        Hess =  @(x) P(dr_adj(d2psi*dr(P(x)))) + H.d2S.d2S(x,H.omega,H.m) + shift*x;
        % Run indirect registration.
        % [v,his] = NPIRpara.scheme(objective_new,v,'yStop',vStop,NO{:});
        dy = pcg(Hess,-dJ',10^-2,250,Preconditioner);
        descent =   dJ * dy;
        if descent > 0
            warning('no descent direction, switch to -dy!')
            dy      = -dy;
        end
        [t,v,LSiter] = Armijo(objective_new,v,dy,J,dJ,...
            'LSMaxIter',10,'LSreduction',1e-2);
    end
    vc=v;
 end
elapsed=toc(ticId);

% Transform template and reshape.
yc = getTrafoFromInstationaryVelocityRK4(vc, getNodalGrid(omega,ML{maxLevel}.m),...
    'omega', omegaV, 'm', ML{maxLevel}.m, 'nt', nt, 'tspan', [1, 0], 'N', N);
rec = imgModel(T, omega, center(yc, ML{maxLevel}.m));
rec = reshape(rec, ML{maxLevel}.m);
z=reshape(z,m);

%save results
[resfile,resfile1,resfile2,resfile3,resfile4,paramfile] = saveres(name, outputfolder, image1,...
    image2,image2, rec,z, dist, reg, objfun, imageModel, N, nt, alpha, theta,...
    noise, elapsed, 1, lambda);
abs_source=abs(z);
imwrite(uint8(255*abs_source/ max(rec(:))), fullfile(outputfolder, sprintf('%s_abs_source_%5f.png', name,norm(abs_source))));
%The deformation error here refers to the error between target and
%deformation part.

