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
function tests = mfThirdOrderSTTest
    tests = functiontests(localfunctions);
end

function resultTest(testCase)

[Sc, dS, d2S] = mfThirdOrderST('para');
verifyEqual(testCase, Sc, 'cell-centered');
verifyEqual(testCase, dS, 1);
verifyEqual(testCase, d2S, @spectralPrecondPCG);

end

function result2dTest(testCase)

% Create zero unknown.
m = [61, 64];
nt = 1;
N = prod(m) * 2 * (nt + 1);
vc = zeros(N, 1);
omega = [0, 1, 0, 1];

alpha = [1, 1, 1];
HessianShift = 1e-2;

[Sc, dS, d2S] = mfThirdOrderST(vc(:), omega, m, 'nt', nt, 'alpha',...
    alpha, 'HessianShift', HessianShift);
verifyEqual(testCase, Sc, 0);
verifyEqual(testCase, dS, zeros(1, N));

% Create unknown.
m = [64, 64];
nt = 1;
N = prod(m) * 2 * (nt + 1);
vc = peaks(64);
vc = repmat(vc(:), 2 * (nt + 1), 1);
omega = [0, 1, 0, 1];

alpha = [1, 1, 1];
HessianShift = 1e-2;

[Sc, dS, d2S] = mfThirdOrderST(vc(:), omega, m, 'nt', nt, 'alpha',...
    alpha, 'HessianShift', HessianShift);
verifyEqual(testCase, Sc > 0, true);
verifyEqual(testCase, size(dS), [1, N]);

end

function MLLDDMMTest(testCase)


% Load images.
image1 = phantom('Modified Shepp-Logan', 64);
image2 = image1;

% Save size of template.
m = size(image1);

% Set bumber of Runge-Kutta steps.
N = 5;

% Define data term.
dist = 'SSD';

% Define regularization term.
reg = 'mfThirdOrderST';

% Define image model.
imageModel = 'splineInterMex';

% Define temporal discretization of the velocity (number of time steps is
% then nt + 1.
nt = 1;

% Set regularization parameters (in order: space, time, L2-norm squared).
% The following parameters can be given as array
alpha = [1, 1, 1];

% Set Hessian shift.
hessianShift = 1e-2;

% Fixed paramters
pad = 0.5;

% Set domain size.
omega = [0, 1, 0, 1];

% Set domain for velocities by padding.
omegaV = omega;
omegaV(1:2:end) = omegaV(1:2:end)-pad;
omegaV(2:2:end) = omega(2:2:end)+pad;

% Initialize models.
imgModel('reset', 'imgModel', imageModel);
regularizer('reset', 'regularizer', reg, 'nt', nt, 'alpha', alpha, 'HessianShift', hessianShift);
trafo('reset', 'trafo', 'affine2D');
distance('reset', 'distance', dist);
viewImage('reset', 'viewImage', 'viewImage2D', 'colormap', flipud(gray(256)));
NPIRpara = optPara('NPIR-GN');
NPIRpara.maxIter = 30;
NPIRpara.scheme = @GaussNewtonLDDMM;

% Create multilevel versions of template.
[ML, minLevel, maxLevel, ~] = getMultilevel({image1, image2}, omega, m, 'fig', false);

% Run algorithm.
mV = @(m) ceil(1*m);
[vc, ~, ~, his] = MLLDDMM(ML, 'minLevel', minLevel, 'maxLevel', maxLevel, 'omegaV', omegaV, 'mV', mV, 'N', N, 'parametric', false, 'NPIRpara', NPIRpara, 'plots', false);

% Transform template and reshape.
yc = getTrafoFromInstationaryVelocityRK4(vc, getNodalGrid(omega,m), 'omega', omegaV, 'm', m, 'nt', nt, 'tspan', [1,0], 'N', N);
Topt = linearInterMex(ML{maxLevel}.T, omega, center(yc, m));
Topt = reshape(Topt, m);

% Output stats.
fprintf('Elapsed time is: %.2f seconds, SSIM=%.2f.\n', his.time, ssim(Topt, image2));

% Verify zero deformation.
verifyEqual(testCase, size(Topt), size(image1));
verifyEqual(testCase, vc, zeros(size(vc)), 'AbsTol', 1e-3);

end