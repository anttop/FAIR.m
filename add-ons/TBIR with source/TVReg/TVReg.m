% This script solves the indirect Tichonov problem 
% 
% min_x 1/2*||K(x)-R||_2^2+ lambda*TV(x)
% 
% for fixed (noisy) data R
% via a primal dual hybrid gradient method with extragradient step
% the forward operator is choosen to be a partial Radon transform.
%

display=true; %display variables
% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'experiments', 'results','TV Reg');
mkdir(outputfolder);

% Load images.
path = fullfile(FAIRpath, 'add-ons', 'experiments', 'data');
file1 = 'deform3.png';
image1 = double(imresize(imread(fullfile(path, file1)), [128,128]));
name = 'phantom';% Name of dataset.
N = size(image1);% Save size of template.


% Set directions for Radon transform.
theta = linspace(0, pi , 40);

% Set up detector size and geometries.
ndet = ceil(hypot(N(1), N(2)));
vol_geom = astra_create_vol_geom(N);
proj_geom = astra_create_proj_geom('parallel', 1.0, ndet, theta);

% Create operator function handles.
K = @(x) radon2d(x, proj_geom, vol_geom);
Kadj = @(x) radon2dadj(x, proj_geom, vol_geom);

% Create cleanup function handle.
cleanup = @() astra_cleanup(proj_geom, vol_geom);
% Set regularisation parameter.
lambda = 10;
noise  = 0.05;
% Create measurements.
R = K(image1);
R = addnoise(R, noise);

X=zeros(N); %initial primal variable
Xbar=X;
A= @(x) my_grad(x,1);         
Aadj= @(x,y) -my_div(x,y,1);  
% initial  dual variable
Yx= zeros(N(1),N(1)); 
Yy= zeros(N(1),N(1));
Yz=zeros(size(R));  

Resolv_F2=@(x,sigma)   perform_soft_thresholding(x, sigma);

%set fix point iteration parameters
tau=1;
sigma1=1/(tau*8) ;
vec=normalize(randn(N));
sigma2=1/(tau*norm(K(vec))^2);

maxIter=10000;
epsilon=0.01;
cond=true; 
count=0;

while cond && count < maxIter
    [temp1, temp2]=A(Xbar);
    temp3=K(Xbar);
    %dual step
    Yx_temp=Yx+sigma1*temp1;
    Yy_temp=Yy+sigma1*temp2;
    Yz_temp=Yz+sigma2*temp3;
    [Yx_new,Yy_new]=Resolv_F1(Yx_temp,Yy_temp,lambda);
    Yz_new=Resolv_F2(Yz_temp,sigma2);
    % primal step
    X_temp=X-tau*Aadj(Yx_new,Yy_new)-tau*Kadj(Yz_new);
    X_new=X_temp;
    Xbar=2*X_new-X;  %extragradient
     
    if norm(X-X_new)/norm(X) < epsilon && norm(Yx-Yx_new)/norm(X) < epsilon ...      
           && norm(Yy-Yy_new)/norm(Yy) < eps && norm(Yz-Yz_new)/norm(Yz) < epsilon
        cond=false;
    end
    X=X_new;
    Yx=Yx_new;
    Yy=Yy_new;
    Yz=Yz_new;
    
     if mod(count,10)==0   && display%display variables
        
       montage({image1/255,R/max(R,[],'all'),X/max(X,[],'all'),Yx/max(abs(Yx),[],'all'),Yy/max(Yy,[],'all'),Yz/max(Yz,[],'all')})
    end
    count=count+1;
end

fprintf('SSIM=%.3f.\n', ssim(X, image1));
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
Rsq = imresize(R, [size(R, 2), size(R, 2)], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino_%2g.png', name,noise)));
imwrite(X/max(max(X(:)), 255) ,  fullfile(outputfolder,sprintf('%s_result_%2g_noise_%2g.png', name,lambda,noise)));


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

function [X,Y]=Resolv_F1(x,y,lamTV)
        nor=(x.^2.+y.^2);
        nor=sqrt(nor);
        if nor==0 ;
            X=x;
            Y=y;
        elseif lamTV==0;
            X=0.*x;
            Y=0.*x;
        else
        temp=max(nor/lamTV,1);
        X=x./temp;
        Y=y./temp;
    end
    
end