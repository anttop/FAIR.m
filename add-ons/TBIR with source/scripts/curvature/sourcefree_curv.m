%This script performes a source-free method to solve a image
%reconstruction problem. The inexact Gauss-Newton method from TBIR is
%adopted. Here, a curvature regulariser is used.

% Flag that activates plotting.
plot = false;

% Set results output folder.
outputfolder = fullfile(FAIRpath, 'add-ons', 'TBIR with source', 'results','curvature', 'source-free method');
mkdir(outputfolder);

% Name of dataset.
name = 'Brain';

% Load images.
path = fullfile(FAIRpath, 'add-ons','TBIR with source',  'data');
file1 = 'brain-T.png';
file2 = 'brain-gausslight0_1.png';
%file3 = 'brain-R.png';
image1 = double(imresize(imread(fullfile(path, file1)), [128, 128]));
image2 = double(imresize(imread(fullfile(path, file2)), [128, 128]));

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
sigma = 0.00;

% Set regularization parameters.
Alpha = [0.001, 0.01, 0.1, 0.2, 0.5,0.7,1];
for i=1:length(Alpha)
   for j=1:length(Alpha)
       alpha=[Alpha(i),Alpha(j)];
       name =append('test',int2str(i),'+', int2str(j));
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
imgModel('reset', 'imgModel', imageModel);
trafo('reset', 'trafo', 'affine2D');
distance('reset', 'distance', dist);
viewImage('reset', 'viewImage', 'viewImage2D', 'colormap', gray(256));
NPIRpara = optPara('NPIR-GN');
NPIRpara.maxIter = 50;
NPIRpara.scheme = @GaussNewtonLDDMM;

% Create multilevel versions of template.
[ML, minLevel, maxLevel, ~] = getMultilevel(image1, omega, m, 'fig', 0);

% Set directions for Radon transform.
theta = linspace(0, 180, 10);
for k=minLevel:maxLevel
    [ML{k}.K, ML{k}.Kadj, ML{k}.cleanup, ML{k}.ndet] = createRadon2d(size(ML{k}.T), theta, 2^(-k + minLevel));
end


% Run algorithm for each setting.
mV = @(m) ceil(1*m);

% Check if measurements exists, otherwise create.
sinogramfile = fullfile(outputfolder, 'Sinograms', sprintf('%s_measurements_%g.mat', name, sigma));
if(exist(sinogramfile, 'file'))
    S = load(sinogramfile);
    ML{maxLevel}.R = S.R;
else
    % Apply operator on finest level to generate synthetic measurements.
    R = ML{maxLevel}.K(image2);
    
    % Add noise to measurements.
    R = addnoise(R, sigma);
    ML{maxLevel}.R = R;

    % Save measurements.
    mkdir(fullfile(outputfolder, 'Sinograms'));
    save(sinogramfile, 'R', 'theta', 'm');
end

% Save template, unknown image, and measurements to results folder.
imwrite(image1 / 255, fullfile(outputfolder, sprintf('%s_source.png', name)));
imwrite(image2 / 255, fullfile(outputfolder, sprintf('%s_target.png', name)));
Rsize = size(ML{maxLevel}.R, 2);
Rsq = imresize(ML{maxLevel}.R, [Rsize, Rsize], 'nearest');
imwrite(Rsq / max(Rsq(:)), fullfile(outputfolder, sprintf('%s_measurements.png', name)));

ML = multilevelRadon2d(ML, maxLevel, minLevel);

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
difference=abs(image2-rec);
nor=norm(difference);
imwrite(difference / 255, fullfile(outputfolder, sprintf('%s_error_result_%5f.png', name,nor)));

% Save result.
[resfile, paramfile] = saveresults(name, outputfolder, image1, image2,...
    ML{maxLevel}.R, rec, dist, reg, objfun, imageModel, N, nt,...
    alpha, theta, sigma, his.time, true);
   end
end
% Free resources.
for k=minLevel:maxLevel
    ML{k}.cleanup();
end
close all;

% Plot result.
if(plot)
    figure;
    colormap gray;
    subplot(1, 4, 1);
    imagesc(image1);
    axis image;
    title('Template image');
    subplot(1, 4, 2);
    imagesc(image2);
    axis image;
    title('Unknown image');
    subplot(1, 4, 3);
    imagesc(ML{maxLevel}.R);
    axis square;
    title('Measurements');
    ylabel('Directions');
    subplot(1, 4, 4);
    imagesc(rec);
    axis image;
    title('Reconstruction');
end
