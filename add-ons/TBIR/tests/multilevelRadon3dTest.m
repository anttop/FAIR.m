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
function tests = multilevelRadon3dTest
    tests = functiontests(localfunctions);
end

function resultTest(testCase)

% Set up test image and projection angles.
D = load('mri.mat', 'D');
data = im2double(squeeze(D.D));
m = size(data);

% Create multilevel versions of image.
[ML, minLevel, maxLevel, ~] = getMultilevel(data, [0, 1, 0, 1, 0, 1], m, 'fig', false);

for k=minLevel:maxLevel
    ML{k}.K = @(x) permute(sum(x, 3), [1, 3, 2]);
    ML{k}.ndet = [ML{k}.m(1), ML{k}.m(2)];
end

% Create measurements on finest level.
ML{maxLevel}.R = ML{maxLevel}.K(data);

% Create multilevel versions of data.
ML = multilevelRadon3d(ML, maxLevel, minLevel);

for k=minLevel:maxLevel-1
    verifyEqual(testCase, size(ML{k}.R, 1), ML{k}.ndet(1));
    verifyEqual(testCase, size(ML{k}.R, 2), 1);
    verifyEqual(testCase, size(ML{k}.R, 3), ML{k}.ndet(2));
end
end

function checkRadon3dDetectorSizeTest(testCase)


% Set up test image and projection angles.
D = load('mri.mat', 'D');
data = im2double(squeeze(D.D));
m = size(data);
theta = 0:10:179;

% Set GPU.
gpuIdx = 0;

% Create multilevel versions of image.
[ML, minLevel, maxLevel, ~] = getMultilevel(data, [0, 1, 0, 1, 0, 1], m, 'fig', false);

% Set up operators for all levels.
for k=minLevel:maxLevel
    [ML{k}.K, ML{k}.Kadj, ML{k}.cleanup, ML{k}.ndet] = createRadon3d(size(ML{k}.T), theta, gpuIdx);
end

try
    % Create measurements on finest level.
    ML{maxLevel}.R = ML{maxLevel}.K(data);

    % Create multilevel versions of data.
    ML = multilevelRadon3d(ML, maxLevel, minLevel);

    % Check if size of downsampled data matches output of Radon operator.
    for k=minLevel:maxLevel-1
        verifyEqual(testCase, size(ML{k}.R), size(ML{k}.K(ML{k}.T)));
    end
catch ME
    warning('%s. Possibly no GPU available or wrong GPU selected.', ME.message);
end

% Free resources.
for k=minLevel:maxLevel
    ML{k}.cleanup();
end
end
