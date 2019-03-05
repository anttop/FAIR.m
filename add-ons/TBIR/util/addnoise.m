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
function y = addnoise(x, sigma)
%ADDNOISE Adds Gaussian white noise of certain level to data.
%
%   y = ADDNOISE(x, sigma) takes an n-dimensional image and adds noise of
%   certain level using imnoise.
%
% Input:
%   x       original data, double array of size m.
%   sigma   noise level (e.g. 0.05 for five percent).
%   
% Output:
%   y       corrupted data, double array of size m.

% Normalise image to [0, 1].
xmax = max(x(:));
xmin = min(x(:));
x = x - xmin;
x = x / xmax;

% Add Gaussian white noise.
y = imnoise(x, 'gaussian', 0, sigma^2 * var(x(:))) * xmax;

end
