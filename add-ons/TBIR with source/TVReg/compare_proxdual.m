%This script compares two implemenation which tackle solve the indirect
%regularisation problem with TV regularisation.
%The first implementation can be found in TVReg.m and the other is used in
%the script TBIR-example 1a). Both utilize a primal dual hybrid gradient
%method.
clear;
close all;
clc;

% Flag that activates plotting.
plot = false;

% toolbox_optim is used for computing TV reconstructions.
addpath(genpath(fullfile(FAIRpath, '../', 'matlab-toolboxes', 'toolbox_optim')));

% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'experiments', 'results', '22.09.2021','TV Reg noise');
mkdir(outputfolder);

% Name of dataset.
name = 'phantom';
% Load images.
path = fullfile(FAIRpath, 'add-ons', 'experiments', 'data');
file1 = 'deform3.png';
file2 = 'deform3.png';
image1 = double(imresize(imread(fullfile(path, file1)), [128,128]));
image2 = double(imresize(imread(fullfile(path, file2)), [128,128]));

% Save size of template.
m = size(image1);

% Set directions for Radon transform.
theta = linspace(0, pi , 10);

% Set up detector size and geometries.
ndet = ceil(hypot(m(1), m(2)));
vol_geom = astra_create_vol_geom(m);
proj_geom = astra_create_proj_geom('parallel', 1.0, ndet, theta);

% Create operator function handles.
K = @(x) radon2d(x, proj_geom, vol_geom);
Kadj = @(x) radon2dadj(x, proj_geom, vol_geom);

% Create cleanup function handle.
cleanup = @() astra_cleanup(proj_geom, vol_geom);

% Create measurements.
R = K(image2);
R = addnoise(R, 0.05);
target=R;
% Set regularisation parameter.
lambda = 0.1;

% Compute TV reconstruction.
K2 = @(x) grad(x);
KS2 = @(x) -div(x);
Amplitude = @(u) sqrt(sum(u.^2, 3));
F2 = @(u) lambda * sum(sum(Amplitude(u)));
Normalize = @(u) u ./ repmat(max(Amplitude(u), 1e-10), [1, 1, 2]);
ProxF2 = @(u, tau) repmat(perform_soft_thresholding(Amplitude(u), lambda * tau), [1 1 2]) .* Normalize(u);
ProxFS = compute_dual_prox(ProxF2);
F = @(x) norm(R - x, 'fro')^2/2;
G = @(x) 0;
ProxF = @(x,tau) (x - tau * R) / (1 + tau / 2);
ProxG = @(x,tau) x;
options.niter = 250;
options.report = @(x) G(x) + F(K(x)) + F2(K2(x));


% Compute TV reconstruction with given template.
[rec2, ~] = perform_dr_pd(zeros(size(Kadj(R))), K, K2, Kadj, KS2, ProxF, ProxFS, ProxG, options);
rec2 = reshape(rec2, m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second implmentation
N=size(image1);
X=zeros(N); %initial primal variable
Xbar=X;

noise=0.05 ;%set noise level

A= @(x) my_grad(x,1);         
Aadj= @(x,y) -my_div(x,y,1);  

Yx= zeros(N(1),N(1)); % initial TV component of dual variable
Yy= zeros(N(1),N(1));
R = @(x) K(x);
Radj =@(x) Kadj(x);
g=target;
Yz=zeros(size(g));   %define HS dual variable
lamTV=lambda;         %set TV reularisation parameter
%Resolv_F1=@(x,sigma)  max(-lamTV,min(x,lamTV));

Resolv_F2=@(x,sigma)  (x-sigma*g) /(1+sigma);
%set fix point iteration parameters
tau=1;
sigma1=1/(tau*8) ;
vec=normalize(randn(N));
sigma2=1/(tau*norm(R(vec))^2);

maxInter=10000;
epsilon=0.001;
cond=true; 
count=0;

while cond & count < maxInter
    [temp1, temp2]=A(Xbar);
    temp3=R(Xbar);
    Yx_temp=Yx+sigma1*temp1;
    Yy_temp=Yy+sigma1*temp2;
    Yz_temp=Yz+sigma2*temp3;
    
    %Yx_new=Resolv_F1(Yx_temp,sigma1);    %dual step
    %Yy_new=Resolv_F1(Yy_temp,sigma1);
    [Yx_new,Yy_new]=Resolv_F1(Yx_temp,Yy_temp,lamTV);
    Yz_new=Resolv_F2(Yz_temp,sigma2);

    
    X_temp=X-tau*Aadj(Yx_new,Yy_new)-tau*Radj(Yz_new);
    X_new=X_temp; % primal step
    Xbar=2*X_new-X;  %extragradient
     
    if norm(X-X_new)/norm(X) < epsilon & norm(Yx-Yx_new)/norm(X) < epsilon ...      
            norm(Yy-Yy_new)/norm(Yy) < eps & norm(Yz-Yz_new)/norm(Yz) < epsilon
        break
    end
    X=X_new;
    Yx=Yx_new;
    Yy=Yy_new;
    Yz=Yz_new;
    count=count+1;
end
rec1=X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print results
fprintf('L2-TV:  SSIM=%.3f.\n', ssim(rec2, image2));
fprintf('my- TV-Reg:, SSIM=%.3f.\n', ssim(rec1, image2));
fprintf('compare results to each other:, SSIM=%.3f.\n', ssim(rec2, rec1));
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(target / 255, fullfile(outputfolder, sprintf('%s_target.png', name)));
Rsq = imresize(g, [size(R, 2), size(R, 2)], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino.png', name)));
imwrite(rec2 / max(max(rec2(:)), 255), fullfile(outputfolder, sprintf('%s_TV-L2_angles_lam_%3f.png', name,lambda)));
imwrite(rec1/max(max(rec1(:)), 255) ,  fullfile(outputfolder,sprintf('%s_myTV_angles_lam_%3f.png', name,lambda)));




% Plot result.
if(plot)
    figure;
    colormap gray;
    subplot(2, 3, 1);
    imagesc(image1, [0, 255]);
    axis image;
    title('Template image');
    subplot(2, 3, 2);
    imagesc(image2, [0, 255]);
    axis image;
    title('Unknown image');
    subplot(2, 3, 3);
    imagesc(R);
    axis square;
    title('Measurements');
    ylabel('Directions');
    subplot(2, 3, 4);
    imagesc(rec1, [0, 255]);
    axis image;
    title('FBP reconstruction');
    subplot(2, 3, 5);
    imagesc(rec2, [0, 255]);
    axis image;
    title('TV reconstruction');
    subplot(2, 3, 6);
    imagesc(rec3, [0, 255]);
    axis image;
    title('TV reconstruction with template');
end

function sino = radon2d(data, proj_geom, vol_geom)
    proj_id = astra_create_projector('linear', proj_geom, vol_geom);
    [sino_id, sino] = astra_create_sino(data, proj_id);
    astra_mex_data2d('delete', sino_id);
    astra_mex_projector('delete', proj_id);
end

function vol = radon2dadj(sino, proj_geom, vol_geom)
    recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    proj_id = astra_create_projector('linear', proj_geom, vol_geom);
    
    cfg = astra_struct('BP');
    cfg.ProjectorId = proj_id;
    cfg.ProjectionDataId = sinogram_id;
    cfg.ReconstructionDataId = recon_id;

    alg_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', alg_id);
    vol = astra_mex_data2d('get', recon_id);

    astra_mex_data2d('delete', sinogram_id, recon_id);
    astra_mex_projector('delete', proj_id);
    astra_mex_algorithm('delete', alg_id);
end

function astra_cleanup(proj_geom, vol_geom)
    astra_mex_data2d('delete', proj_geom, vol_geom);
end

function [x,R] = perform_dr_pd(x, K, K2, KS, KS2, ProxFS, ProxFS2, ProxG, options)

options.null = 0;
report = getoptions(options, 'report', @(x)0);
niter = getoptions(options, 'niter', 100);

if isnumeric(K)
    K = @(x)K*x;
end  
if isnumeric(KS)
    KS = @(x)KS*x;
end      
%%%% Douchglas Rachford parameters %%%%
sigma = getoptions(options, 'sigma', -1);
tau   = getoptions(options, 'tau', -1);
if sigma<0 || tau<0
    [L,~] = compute_operator_norm(@(x) KS(K(x)),randn(size(x)));
    [L2,~] = compute_operator_norm(@(x) KS2(K2(x)),randn(size(x)));
    sigma = .99/(L);
    sigma2 = .99/(L2);
    tau = 1;
end

y = K(x);
y2 = K2(x);
clear R;
lam_k = 1.5;

for i=1:niter
    R(i) = report(x);
    tmp_domain = KS(y);
    p1 = KS2(y2);
    tmp_domain = tmp_domain + p1;
    tmp_domain = x -tau/2*tmp_domain;
    p1 = ProxG(tmp_domain,tau);
    w1 = 2*p1 -x;
    tmp = y + (sigma/2) * K(w1);
    p2 = ProxFS(tmp, sigma);
    w2 = 2*p2 -y;
    tmp = y2 + (sigma2/2) * K2(w1);
    p22 = ProxFS2(tmp, sigma2);
    w22 = 2*p22 -y2;
    tmp_domain = KS(w2);
    z1= KS2(w22);
    tmp_domain = tmp_domain + z1;
    z1 = w1 - (tau/2)*tmp_domain;
    x = x + lam_k*(z1-p1);
    tmp_domain = 2*z1 - w1;
    z2 = w2 + (sigma/2)*K(tmp_domain);
    y = y + lam_k*(z2 -p2);
    z22 = w22 + (sigma2/2)*K2(tmp_domain);
    y2 = y2 + lam_k*(z22 -p22);
end
x = p1;
end

function [X,Y]=Resolv_F1(x,y,lamTV)
        nor=(x.^2.+y.^2);
        nor=sqrt(nor);
        if nor==0 ;
            X=x;
            Y=y;
        else
        temp=max(nor/lamTV,1);
        X=x./temp;
        Y=y./temp;
    end
    
end
    
