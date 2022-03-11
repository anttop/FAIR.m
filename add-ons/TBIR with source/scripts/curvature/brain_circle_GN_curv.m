%This script is a modified version of the script 'brain_circle_curv.m'.
%Instead of the iPALM scheme, we used a method which combines the PALM
%algorithm with the inexact Gauss-Newton scheme found in 'LDDMM.m'.
%Again, starting from initial values (v,z), in every iteration we perform
%i)   A inexact Gauss-Newton step on v for fixed z
%)ii) and update z via a first oder step (PALM).
% Inertia and postprocessing are dropped.
loadData = false;
plot=false;
name1='brain_GN';
% Load images.
path = fullfile(FAIRpath, 'add-ons', 'TBIR with source', 'data');
file1 = 'brain-T.png';
file2 = 'brain-circle.png';
file3 = 'brain-R.png';
image1 = double(imresize(imread(fullfile(path, file1)), [128, 128]));
image2 = double(imresize(imread(fullfile(path, file2)), [128, 128]));
image3 = double(imresize(imread(fullfile(path, file3)), [128, 128]));
m = size(image1);

% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'TBIR with source', 'results','curvature','brain circle', 'Newton step');
mkdir(outputfolder);
% Save template, unknown image, and measurements to results folder
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name1)));
imwrite(image2 / 255, fullfile(outputfolder, sprintf('%s_target.png', name1)));

% Set parameters and models
N = 5; % Runge-Kutta steps
nt = 1; % Temporal discretization of velocity
noise = 0.05; % Noise
alpha = [1,1]; % Regularization velocity.
lambda = 0.2; % Regularization source
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
shift = 1e-3;
regularizer('reset', 'regularizer', reg, 'nt', nt,...
    'alpha', alpha, 'HessianShift', shift);
objfun1 = str2func('LDDMMobjFctn');
NPIRpara = optPara('NPIR-GN');
NPIRpara.maxIter = 2;
NPIRpara.scheme = @GaussNewtonLDDMM;
NPIRpara.Plots = @FAIRplots_op;
NO = FAIRcell2struct(NPIRpara);
    
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


Rsize = size(ML{maxLevel}.R, 2);
Rsq = imresize(ML{maxLevel}.R, [Rsize, Rsize], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino_%.2f.png', name, noise)));

%PALM
ProxG = @(x,sigma) proxTV_toolbox(x,sigma*lambda); %second variable's 
maxIter=100;
tol=10^-4;
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
    
    % Initialize run
    v=v0; % Starting guess
    vi=v;
    zi=z;
    J=0;
    iter=1;
    Knorm=operator_norm(K,Kadj,rand(m(1)*m(2),1));
    sigma=Knorm^2; %step size second variable
    tau=sigma; % Equal parameters for step sizes
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
        
        % Newton update v with line search
        Rv=Rc-K(z);
        objective_new = @(vc) objfun1(T,Rv,omega,m,vRef,xc,omegaV,mV(m),N,vc,K,Kadj);
        [J,~,dJ,H] = objective_new(v);
        Amplitude = @(u) sqrt(sum(u.^2, 3));
        J + lambda * sum(sum(Amplitude(grad(z)))) %monitor objective value 
        P       = H.d2D.P; % grid-to-grid operator
        dr      = H.d2D.dr;
        dr_adj  = H.d2D.dr_adj;
        d2psi   = H.d2D.d2psi;
        Hess =  @(x) P(dr_adj(d2psi*dr(P(x)))) + H.d2S.d2S(x,H.omega,H.m) + shift*x;
        % Run indirect registration.
        % [v,his] = NPIRpara.scheme(objective_new,v,'yStop',vStop,NO{:});
        dy = pcg(Hess,-dJ',10^-2,250,Preconditioner);
        descent =   dJ * dy;
        if descent > 0,
            warning('no descent direction, switch to -dy!')
            dy      = -dy;
        end
        [t,v,LSiter] = Armijo(objective_new,v,dy,J,dJ,...
        'LSMaxIter',10,'LSreduction',1e-2);
        
        % Update z (second variable)
        yc = getTrafoFromInstationaryVelocityRK4(v, getNodalGrid(omega,m),...
            'omega', omegaV, 'm', m, 'nt', nt, 'tspan', [1 0], 'N', N);
        [rec,~] = imgModel(T, omega, center(yc, m));
        Rz=Rc-K(reshape(rec,m));
         for ii=1:10
             dD= reshape(Kadj((K(z)-Rz)),m);
            z_temp=z-dD./sigma;
            z=ProxG(z_temp,sigma);
         end
        
        % Create plot
        A=reshape(rec, m);
        AA=A+z;
        %'maxdef', max(rec), 'max z' ,max(z(:)), 'max rec' , max(AA(:))
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
    vc = v;
end
elapsed=toc(ticId);

% Transform template and reshape.
yc = getTrafoFromInstationaryVelocityRK4(vc, getNodalGrid(omega,ML{maxLevel}.m),...
    'omega', omegaV, 'm', ML{maxLevel}.m, 'nt', nt, 'tspan', [1, 0], 'N', N);
rec = imgModel(T, omega, center(yc, ML{maxLevel}.m));
rec = reshape(rec, ML{maxLevel}.m);
z=reshape(z,ML{maxLevel}.m);
%save results
[resfile,resfile1,resfile2,resfile3,resfile4,paramfile] = saveres(name, outputfolder, image1,...
    image2,image3, rec,z, dist, reg, objfun, imageModel, N, nt, alpha, theta,...
    noise, elapsed, 0, lambda);