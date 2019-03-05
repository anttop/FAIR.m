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
function [D, r, dD, dr, dradj, d2psi] = SSD_op(Tc, Rc, omega, m, K, Kadj, varargin)
%SSD_OP Computes sum-of-squares objective including an operator.
%
%   [D, r, dD, dr, dr_adj, d2psi] = SSD_OP(Tc, Rc, omega, m, K, Kadj, varargin)
%   takes a deformed image Tc, measurements Rc, the size of the measurement
%   domain omega, the number of sampling points m of the measurement domain,
%   and arguments varargin, and returns the value of the evaluated
%   sum-of-squared distances data term and derivatives.
%
%   D = SSD(K(), Rc) = hd * psi(r(Tc))
%   
%   with r(Tc) = K(Tc) - Rc and psi(r) = 0.5 * r' * W * r.
%
%   The first derivative of r(Tc) is dr = W * r,
%   the second derivative of psi(r) is d2psi = hd * weights,
%   and hd the cell size.
%
% Input:
%   Tc          array of deformed template and of size n.
%   Rc          array of measurements and of size m.
%   omega       vector denoting the size of the
%               measurement geometry (e.g. [0, 1] for Radon transform of a
%               2D image.
%   m           is a vector of dimension n - 1 denoting the number of
%               sample points of the measurement geometry.
%   K           linear operator mapping from R^prod(n) to R^prod(m).
%   Kadj        adjoint of linear operator K.
%   varargin    optional name-value pairs such as 'doDerivative' (boolean)
%               or a scalar weights or a (diagonal) weight matrix of
%               size [prod(m), prod(m)].
% Output:
%   D           non-negative scalar (provided W is non-negative).
%   r           array of size m.
%   dD          array of size n.
%   dr          array of size m.
%   dr_adj      array of size 
%   d2psi       array of same size as weights if set, otherwise a scalar.
%
dD = [];
dr = [];
dradj = [];
doDerivative = false;
weights = 1.0;

% Parse parameters.
for k=1:2:length(varargin)
  eval([varargin{k}, '=varargin{', int2str(k+1), '};']);
end

% Compute weights matrix.
hd = prod((omega(2:2:end)-omega(1:2:end)) ./ m);
d2psi = hd * weights;

% Compute residual.
r = K(Tc) - Rc;

% Compute SSD with operator.
D = 0.5 * (r(:)' * d2psi * r(:));

if(doDerivative)
    dr = K;
    dradj = Kadj;
    dD = Kadj(d2psi * r)';
else
end
