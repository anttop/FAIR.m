%% This script solves the proposed problem with a curvature regulariser instead of
% third order regularisation.
__________________________________________________________________________
loadData = false;
plot=false;
name='hands';
% Load images.
path = fullfile(FAIRpath, 'kernel',  'data');
file1 = 'hands-R.jpg';
file2 = 'Hands-T.jpg';
file3 = 'Hands-T.jpg';
image1 = double(imresize(imread(fullfile(path, file1)), [128, 128]));
image2 = double(imresize(imread(fullfile(path, file2)), [128, 128]));
image3 = double(imresize(imread(fullfile(path, file3)), [128, 128]));
m = size(image1);
% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'TBIR with source', 'results','curvature','hands');
mkdir(outputfolder);
% Save template, unknown image, and measurements to results folder
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(image2 / 255, fullfile(outputfolder, sprintf('%s_target.png', name)));


% Set parameters and models
N = 5; % Runge-Kutta steps
nt = 1; % Temporal discretization of velocity
noise = 0.05; % Noise
alpha = [6,6]; % Regularization velocity.
lambda = 1; % Regularization source
dist = 'SSD_op';
reg = 'mfCurvatureST';
objfun = 'PALMobjFctn';
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

% Create multilevel versions of data and operators
[ML, ~, maxLevel, ~] = getMultilevel(image1, omega, m, 'fig', 0);
minLevel=4;

% Set directions for Radon transform.
theta = linspace(0, 180, 10);
for k=minLevel:maxLevel
    [ML{k}.K, ML{k}.Kadj, ML{k}.cleanup, ML{k}.ndet] = createRadon2d(size(ML{k}.T), theta, 2^(-k + minLevel));
end

% Load or create measurements
sinogramfile = fullfile(outputfolder, 'Sinograms', sprintf('%s_sino_%g_%g.mat', name, noise, length(theta)));
if loadData && (exist(sinogramfile, 'file'))
    S = load(sinogramfile);
    ML{maxLevel}.R = S.R;
else
    % Create and save measurements
    R = ML{maxLevel}.K(image2);
    R = addnoise(R, noise);
    ML{maxLevel}.R = R;
    mkdir(fullfile(outputfolder, 'Sinograms'));
    save(sinogramfile, 'R', 'theta', 'm');
end
% Create multilevel versions of measurements.
ML = multilevelRadon2d(ML, maxLevel, minLevel);

%save sinogram
Rsize = size(ML{maxLevel}.R, 2);
Rsq = imresize(ML{maxLevel}.R, [Rsize, Rsize], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino_%.2f.png', name, noise)));

%PALM 
ProxG = @(x,sigma) proxTV_toolbox(x,sigma*lambda); %proximal mapping for second variable
maxIter=100;
tol=10^-4;
vRef=[];
eta=3; % Backtracking parameter
malLevel = 7;

ticId = tic; 
for level=minLevel:maxLevel
    % Update m,grid,data and coefficients
    m    = ML{level}.m;
    xc   = getCellCenteredGrid(omega,m);
   [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
    Rc = ML{level}.R;
    K=ML{level}.K;
    Kadj=ML{level}.Kadj;
    trafo([],xc);
    
    % compute starting guess y0 for this level and stopping value yStop
    if level == minLevel
        if isempty(vRef)
            vRef = getVelocityStartingGuess(omegaV,mV(m),nt);
        end
        v0 = vRef;  % the best known so far
        z = zeros(m);
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
        z = imresize(z,m,'nearest');
    end
    vStop = 0*getVelocityStartingGuess(omegaV,mV(m),nt);
    
    % Initialize iPALM run
    v=v0; % Starting guess
    vi=v;
    zi=z;
    J=0;
    iter=1;
    Knorm=operator_norm(K,Kadj,rand(m(1)*m(2),1));
    sigma=Knorm^2;
    tau=sigma; % Equal parameters for stepsizes
    
    %proximal mapping for first variable
    [~,~,ScdDev] = regularizer(v0,omegaV,m,'doDerivative',1,'tspan',[1 0]);
    ProxF = @(x)x + ScdDev.d2S(x,omega,m)/tau;
    Preconditioner = @(x) solveSpectral(omega,m,[1 0],nt,1,1/tau,x);
    STOP(1) = false;
    STOP(2) = false;
    count=0;
    
    while ~any(STOP)
        % Keep record
        v_old = v;
        z_old = z;
        J_old = J;
        
       % Update v with line search
        Rv=Rc-K(z);
        NPIRfctn = @(vc) PALMobjFctn(T,Rv,omega,m,vRef,xc,omegaV,mV(m),N,vc,K,Kadj);
        J=NPIRfctn(v);
        [~,~,dJ,~]=NPIRfctn(vi);
        
        %determine stepsize tau via backtracking scheme
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
        Amplitude = @(u) sqrt(sum(u.^2, 3));
        J1+lambda * sum(sum(Amplitude(grad(z)))); %monitor function value
        v = v_temp;
        vi = v + (iter-1)/(iter+2).*(v-v_old);
        
        % Update z without linesearch
        yc = getTrafoFromInstationaryVelocityRK4(v, getNodalGrid(omega,m),...
            'omega', omegaV, 'm', m, 'nt', nt, 'tspan', [1 0], 'N', N);
        [rec,~] = imgModel(T, omega, center(yc, m));
         Rz=Rc-K(reshape(rec,m));
         dD= reshape(Kadj((K(zi)-Rz)),m);
         z_temp=zi-dD./sigma;
         z=ProxG(z_temp,sigma);
         zi = z + (iter-1)/(iter+2).*(z-z_old);
        
        % Create plot
        A=reshape(rec, m); %deformation part
        AA=A+z;            %complete reconstruction
       % 'maxdef', max(rec), 'max z' ,max(z(:)), 'max rec' , max(AA(:));
        if plot && mod(iter,10)==0 
          %display images
          montage({rescale(A),z/(max(abs(z(:)))+0.001), rescale(AA)});
        end
        % Set variables for next iteration
        iter = iter + 1;
        STOP(1) = iter > maxIter;
        STOP(2) = (norm(v-v_old) < tol*norm(v)) && (norm(z-z_old) < tol*norm(z));
     %   STOP(3) = (abs(J-J_old) < tol*abs(J));
    end
    
    
    %postprocessing: inexact Gauss-Newton on variable v
    % Define objective.
    objfun1 = 'LDDMMobjFctn';
     % Initialize models.
    imgModel('reset', 'imgModel', imageModel);
    trafo('reset', 'trafo', 'affine2D');
    distance('reset', 'distance', dist);
    viewImage('reset', 'viewImage', 'viewImage2D', 'colormap', gray(256));
    NPIRpara = optPara('NPIR-GN');
    NPIRpara.maxIter = 50;
    NPIRpara.scheme = @GaussNewtonLDDMM;
    hessianShift = 1e-2;
   
    ML{level}.R=Rc-K(z);
     % Run indirect registration.
    regularizer('reset', 'regularizer', reg, 'nt', nt,...
        'alpha', alpha, 'HessianShift', hessianShift);
    [vc, ~, ~, his] = MLLDDMM(ML, 'operator', true, 'minLevel',...
        level, 'maxLevel', level, 'omegaV', omegaV, 'mV', mV,...
        'N', N, 'parametric', false, 'NPIRpara', NPIRpara,...
        'NPIRobj', str2func(objfun1), 'plots', plot,'vRef', v);
end

elapsed=toc(ticId);
% Transform template and reshape.
yc = getTrafoFromInstationaryVelocityRK4(vc, getNodalGrid(omega,ML{maxLevel}.m),...
    'omega', omegaV, 'm', ML{maxLevel}.m, 'nt', nt, 'tspan', [1, 0], 'N', N);
rec = imgModel(T, omega, center(yc, ML{maxLevel}.m));
rec = reshape(rec, ML{maxLevel}.m);

%save results
[resfile,resfile1,resfile2,resfile3,resfile4,paramfile] = saveres(name, outputfolder, image1,...
    image2,image3, rec,z, dist, reg, objfun, imageModel, N, nt, alpha, theta,...
    noise, elapsed, 0, lambda);
