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
% This script creates the results shown in Figure 9.
% This script requires data from https://www.fips.fi/dataset.php
%
% Dataset available at https://doi.org/10.5281/zenodo.1183532 (Version 1)
%
% Download and place the following files in the folder 'data':
%
%   - DataDynamic_128x60.mat
%
clear;
close all;
clc;

% Flag that activates plotting.
plot = false;

% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'TBIR', 'results', 'ex_8');
mkdir(outputfolder);

% Name of dataset.
name = 'dynamic';

% Load measurements.
D = load('DataDynamic_128x60.mat');
sinogram = D.sinogram;

% Extract measurements for one time instant and set directions.
t = 1;
slice = sinogram(:, (t-1)*60+1:t*60);

% Reconstruct image using FBP as template.
D = 2240 * 54 / 12;
FSS = 54 / 63;
G = 2240 * (63 - 54)/12;
image1 = ifanbeam(slice, D,...
    'FanRotationIncrement', 1,...
    'FanSensorGeometry', 'line',...
    'FanSensorSpacing', FSS,...
    'OutputSize', 128);

% Extract measurements for one time instant and set directions.
t = 4;
slice = sinogram(:, (t-1)*60+1:t*60);

% Reconstruct image using FBP as ground truth.
D = 2240 * 54 / 12;
FSS = 54 / 63;
image2 = ifanbeam(slice, D,...
    'FanRotationIncrement', 1,...
    'FanSensorGeometry', 'line',...
    'FanSensorSpacing', FSS,...
    'OutputSize', 128);

% Set directions for Radon transform.
theta = 1+(t-1):6:360+(t-1);

% Save size of template.
m = size(image1);

% Set bumber of Runge-Kutta steps.
N = 5;

% Define data term.
dist = 'NCC_op';

% Define regularization term.
reg = 'mfThirdOrderST';

% Define image model.
imageModel = 'splineInterMex';

% Define temporal discretization of the velocity (number of time steps is
% then nt + 1.
nt = 1;

% Define noise level.
sigma = 0;

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
distance('reset', 'distance', dist);
imgModel('reset', 'imgModel', imageModel);
trafo('reset', 'trafo', 'affine2D');
viewImage('reset', 'viewImage', 'viewImage2D', 'colormap', gray(256));
NPIRpara = optPara('NPIR-GN');
NPIRpara.maxIter = 50;
NPIRpara.scheme = @GaussNewtonLDDMM;

% Create multilevel versions of template.
[ML, ~, maxLevel, ~] = getMultilevel(image1, omega, m, 'fig', 0);

% Set starting level.
minLevel = maxLevel - 3;

% Set up operators for all levels.
ndet = size(slice, 1);
for k=minLevel:maxLevel
    [ML{k}.K, ML{k}.Kadj, ML{k}.cleanup, ML{k}.ndet] = createRadon2d_fan(size(ML{k}.T), theta, ceil(ndet * 2^(k - maxLevel)), 2^(-k + minLevel), D, G);
end

% Run algorithm for each setting.
mV = @(m) ceil(1*m);

% Set measurements.
ML{maxLevel}.R = permute(slice, [2, 1]);

% Save template, unknown image, and measurements to results folder.
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(image2 / 255, fullfile(outputfolder, sprintf('%s_target.png', name)));
Rsize = size(ML{maxLevel}.R, 2);
Rsq = imresize(ML{maxLevel}.R, [Rsize, Rsize], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_sino.png', name)));

% Create multilevel versions of measurements.
ML = multilevelRadon2d_imresize(ML, maxLevel, minLevel);

% Define objective.
objfun = 'LDDMMobjFctn';

% Set regularization parameters.
alpha = [2.5e-2, 1, 1e-6];

% Run indirect registration.
regularizer('reset', 'regularizer', reg, 'nt', nt,...
    'alpha', alpha, 'HessianShift', hessianShift);
[vc, ~, ~, his] = MLLDDMM(ML, 'operator', true, 'minLevel',...
    minLevel, 'maxLevel', maxLevel, 'omegaV', omegaV, 'mV', mV,...
    'N', N, 'parametric', false, 'NPIRpara', NPIRpara,...
    'NPIRobj', str2func(objfun), 'plots', plot);

% Transform template and reshape.
yc = getTrafoFromInstationaryVelocityRK4(vc, getNodalGrid(omega,m),...
    'omega', omegaV, 'm', m, 'nt', nt, 'tspan', [1, 0], 'N', N);
rec = linearInterMex(ML{maxLevel}.T, omega, center(yc, m));
rec = reshape(rec, m);

% Output stats.
fprintf('Elapsed time is: %.2f seconds, SSIM=%.3f.\n', his.time, ssim(rec, image2));

% Save result.
[resfile, paramfile] = saveresults(name, outputfolder, image1, image2,...
    ML{maxLevel}.R, rec, dist, reg, objfun, imageModel, N, nt,...
    alpha, theta, sigma, his.time);

% Free resources.
for k=minLevel:maxLevel
    ML{k}.cleanup();
end
close all;

% Plot result.
if(plot)
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
    imagesc(ML{maxLevel}.R);
    axis square;
    title('Measurements');
    ylabel('Directions');
    subplot(2, 3, 4);
    imagesc(rec);
    axis image;
    title('NCC, transport equation');
end
