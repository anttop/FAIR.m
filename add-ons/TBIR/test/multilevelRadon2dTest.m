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
function tests = multilevelRadon2dTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function oddSizeTest(testCase)

% Create multilevel versions of image.
img = phantom('Modified Shepp-Logan', 99);
[ML, minLevel, maxLevel, ~] = getMultilevel(img, [0, 1, 0, 1], size(img), 'fig', 0);

% Set data.
ML{maxLevel}.R = img;

% Create multilevel versions of data.
ML = multilevelRadon2d(ML, maxLevel, minLevel);

for k=minLevel:maxLevel-1
    verifyEqual(testCase, size(ML{k}.R, 2), ceil(size(ML{k+1}.R, 2) / 2));
end
end

function evenSizeTest(testCase)

% Create multilevel versions of image.
img = phantom('Modified Shepp-Logan', 100);
[ML, minLevel, maxLevel, ~] = getMultilevel(img, [0, 1, 0, 1], size(img), 'fig', 0);

% Set data.
ML{maxLevel}.R = img;

% Create multilevel versions of data.
ML = multilevelRadon2d(ML, maxLevel, minLevel);

for k=minLevel:maxLevel-1
    verifyEqual(testCase, size(ML{k}.R, 2), ceil(size(ML{k+1}.R, 2) / 2));
end
end
