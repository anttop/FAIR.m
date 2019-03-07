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
function tests = createRadon3dTest
    tests = functiontests(localfunctions);
end

function operatorTest(testCase)

% Set up test image and projection angles.
D = load('mri.mat', 'D');
data = im2double(squeeze(D.D));
m = size(data);
theta = 0:10:179;

% Set GPU.
gpuIdx = 0;

% Create operators.
[K, Kadj, cleanup, ndet] = createRadon3d(m, theta, gpuIdx);
verifyEqual(testCase, ndet, [m(3), ceil(sqrt(sum(m(1:2).^2)))]);

try
    % Apply operator and check size.
    y = K(data(:));
    verifyEqual(testCase, size(y), [ceil(sqrt(sum(m(1:2).^2))), length(theta), m(3)]);

    % Apply adjoint and check size.
    xrecon = Kadj(y);
    verifyEqual(testCase, size(xrecon), size(data(:)));
catch ME
    warning('%s. Possibly no GPU available or wrong GPU selected.', ME.message);
end

% Free resources.
cleanup();

end

function adjointnessTest(testCase)

% Set up random test data.
data = rand(100, 100, 40);
m = size(data);
theta = 0:10:179;
n = [ceil(sqrt(sum(m(1:2).^2))), length(theta), m(3)];
meas = rand(n);

% Set GPU.
gpuIdx = 0;

% Create operators.
[K, Kadj, cleanup, ndet] = createRadon3d(m, theta, gpuIdx);
verifyEqual(testCase, ndet, [m(3), ceil(sqrt(sum(m(1:2).^2)))]);

try
    % Verify adjointness.
    Kx = K(data(:));
    Kadjy = Kadj(meas);
    verifyEqual(testCase, abs(Kx(:)'*meas(:) - data(:)'*Kadjy), 0, 'AbsTol', 1e-3);
catch ME
    warning('%s. Possibly no GPU available or wrong GPU selected.', ME.message);
end
    
% Free resources.
cleanup();

end
