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
function img = createdisk(m, c, r)
%CREATEDISK Creates a 2D disk image.
%
%   img = CREATEDISK(m, n, c, r) takes an image size m, a centre c, and a
%   radius r, and returns an image of size m with a disk of radius r placed
%   around c.
%
% Input:
%   m       a vector of length two.
%   c       a vector of length two.
%   r       a scalar > 0.
%   
% Output:
%   img     an array of size m.

% Create image.
img = zeros(m);
[X, Y] = ndgrid(1:m(1), 1:m(2));
idx = hypot(X - c(1), Y - c(2)) <= r;
img(idx) = 1;

end
