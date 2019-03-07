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
function ML = multilevelRadon3d(ML, maxLevel, minLevel)
%MULTILEVELRADON3D Creates multlevel versions of 3D Radon transform
%measurements.
%
%   ML = MULTILEVELRADON3D(ML, maxLevel, minLevel) takes a multilevel cell
%   array generated by getMultilevel and creates appropriate multilevel
%   versions of the 2D measurements generated by the 3D Radon transform by
%   simple downsampling.
%
% Input:
%   ML          cell array generated by getMultilevel.
%   maxLevel    the level such that ML{maxLevel}.R are the original
%               measurements.
%   minLevel    the coarsest level at which measurements are generated.
%
% Output:
%   ML          cell array with measurements ML{k}.R generated for all
%               levels between minLevel and maxLevel.
%
% Note that ML{k}.ndet must be set for each level k so that the data can be
% appropriately downsampled.
%
% Note that ML{maxLevel}.R stays untouched.

for k=1:maxLevel - minLevel
    % Current level.
    level = maxLevel-k;
    
    % Select data from higher level.
    data = ML{level + 1}.R;
    
    % Resize data according to number of detectors.
    [~, n, ~] = size(data);
    newsize = [ML{level}.ndet(1), n, ML{level}.ndet(2)];
    ML{level}.R = zeros(newsize);
    for l=1:n
        ML{level}.R(:, l, :) = imresize(squeeze(data(:, l, :)), [ML{level}.ndet(1), ML{level}.ndet(2)], 'bilinear', 'Antialiasing', false) / 2;
    end
end
end
