% Copyright 2019 Lukas F. Lang and Sebastian Neumayer
%
% This file is part of TBIR.
%
%    TBIR is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    TBIR is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with TBIR.  If not, see <http://www.gnu.org/licenses/>.
%
% This script requires the toolbox_optim:
%
% matlab-toolboxes: Matlab toolboxes from www.numerical-tours.com.
% GitHub: https://github.com/gpeyre/matlab-toolboxes
% URL: http://www.numerical-tours.com/
% Version used: 0cd622c
% 
% Clone from github e.g. to some directory by:
%
% >> git clone https://github.com/gpeyre/matlab-toolboxes.git
% 
% and set the path below accordingly. In order to work with the
% abovementioned version type:
% 
% >> cd matlab-toolboxes
% >> git checkout 0cd622c
%
% This script creates the results shown in Figure 1.
clear;
close all;
clc;

% toolbox_optim is used for computing TV reconstructions.
addpath(genpath(fullfile(FAIRpath, '../', 'matlab-toolboxes', 'toolbox_optim')));

% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'TBIR', 'results', 'Ex_1');
mkdir(outputfolder);

% Name of dataset.
name = 'Brain';

% Load images.
path = fullfile(FAIRpath, 'add-ons', 'TBIR', 'data');
file1 = 'brain-T.png';
file2 = 'brain-R.png';
image1 = double(imread(fullfile(path, file1)));
image2 = double(imread(fullfile(path, file2)));

% Save size of template.
m = size(image1);

% Set directions for Radon transform.
theta = linspace2(0, pi / 3, 6);

% Set up detector size and geometries.
ndet = ceil(hypot(m(1), m(2)));
vol_geom = astra_create_vol_geom(m);
proj_geom = astra_create_proj_geom('parallel', 1.0, ndet, theta);

% Create operator function handles.
K = @(x) radon2d(x, proj_geom, vol_geom);
Kadj = @(x) radon2dadj(x, proj_geom, vol_geom);

% Create cleanup function handle.
cleanup = @() astra_cleanup(proj_geom, vol_geom);

% Create measurements and add noise.
R = K(image2);

% Set regularisation parameter.
lambda = 1;

% Compute filtered backprojection.
ticId = tic;
rec1 = iradon(radon(image2, 180 * theta / pi), 180 * theta / pi, 'linear', 'Ram-Lak', 1, m(1));
elapsed1 = toc(ticId);

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
ticId = tic;
[rec2, ~] = perform_dr_pd(zeros(size(Kadj(R))), K, K2, Kadj, KS2, ProxF, ProxFS, ProxG, options);
elapsed2 = toc(ticId);
rec2 = reshape(rec2, m);

template = grad(image1);
F2 = @(u) lambda * sum(sum(Amplitude(u - template)));
Normalize = @(u) u ./ repmat(max(Amplitude(u), 1e-10), [1, 1, 2]);
ProxF2 = @(u, tau) repmat(perform_soft_thresholding(Amplitude(u - template), lambda * tau), [1, 1, 2]) .* Normalize(u - template) + template;
ProxFS = compute_dual_prox(ProxF2);
ticId = tic;
[rec3, ~] = perform_dr_pd(zeros(size(Kadj(R))), K, K2, Kadj, KS2, ProxF, ProxFS, ProxG, options);
elapsed3 = toc(ticId);
rec3 = reshape(rec3, m);

% Free ASTRA resources.
cleanup();

% Output and save results.
fprintf('FBP: elapsed time is: %.2f seconds, SSIM=%.3f.\n', elapsed1, ssim(rec1, image2));
fprintf('L2-TV: elapsed time is: %.2f seconds, SSIM=%.3f.\n', elapsed2, ssim(rec2, image2));
fprintf('L2-TV-template: elapsed time is: %.2f seconds, SSIM=%.3f.\n', elapsed3, ssim(rec3, image2));
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(image2 / 255, fullfile(outputfolder, sprintf('%s_target.png', name)));
Rsq = imresize(R, [size(R, 2), size(R, 2)], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino.png', name)));
imwrite(rec1 / max(rec1(:)), fullfile(outputfolder, sprintf('%s_result_FBP.png', name)));
imwrite(rec2 / max(rec2(:)), fullfile(outputfolder, sprintf('%s_result_L2TV.png', name)));
imwrite(rec3 / max(rec3(:)), fullfile(outputfolder, sprintf('%s_result_L2TV_with_template.png', name)));

% Plot result.
figure;
colormap gray;
subplot(2, 3, 1);
imagesc(image1);
axis image;
title('Template image');
subplot(2, 3, 2);
imagesc(image2);
axis image;
title('Unknown image');
subplot(2, 3, 3);
imagesc(R);
axis square;
title('Measurements');
ylabel('Directions');
subplot(2, 3, 4);
imagesc(rec1);
axis image;
title('FBP reconstruction');
subplot(2, 3, 5);
imagesc(rec2);
axis image;
title('TV reconstruction');
subplot(2, 3, 6);
imagesc(rec3);
axis image;
title('TV reconstruction with template');

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
