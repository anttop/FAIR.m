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
function [K, Kadj, cleanup, ndet] = createIdentity2d(m, scale)
%CREATEIDENTITY2D Creates function handles for the identity operator.
%
%   [K, Kadj, cleanup, ndet] = CREATEIDENTITY2D(m, theta, scale) takes size
%   of 2D geometry m, and a scaling factor, and returns function handles.
%
%   Note that both K and Kadj accept and return vectors!
%
% Input:
%   m       size of 2D geometry, vector of length two.
%   scale   a scalar > 0.
%   
% Output:
%   K       function handle, K maps a vector of length prod(m) to a matrix
%           of size [length(theta), ndet].
%   Kadj    function handle, Kadj maps a matrix of size
%           [length(theta), ndet] to a vector of length prod(m).

% Create operator function handles.
K = @(x) scale * reshape(x, m);
Kadj = @(x) reshape(scale * x, [], 1);
cleanup = @(x) 1;
ndet = m;

end
