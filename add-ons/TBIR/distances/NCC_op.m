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
function [D, r, dD, dr, dradj, d2psi] = NCC_op(Tc, Rc, omega, m, K, Kadj, varargin)
%NCC_OP Computes normalised cross-correlation including an operator.
%
%   [D, r, dD, dr, dr_adj, d2psi] = NCC_OP(Tc, Rc, omega, m, K, Kadj, varargin)
%   takes a deformed image Tc, measurements Rc, the size of the measurement
%   domain omega, the number of sampling points m of the measurement domain,
%   and arguments varargin, and returns the value of the evaluated
%   NCC distance data term and derivatives.
%
%   D = NCC(K(), Rc) = hd * psi(r(Tc))
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
%   varargin    optional name-value pairs such as 'doDerivative' (boolean).
%
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
d2psi = [];
doDerivative = false;

% Parse parameters.
for k=1:2:length(varargin)
  eval([varargin{k}, '=varargin{', int2str(k+1), '};']);
end

% Scaling factor.
s = 1e6;

% Compute first argument.
r = K(Tc);

% Compute second argument.
y = Rc / norm(Rc, 'fro');

% Compute squared norm of r.
snr = r(:)' * r(:);

% Compute NCC with operator.
ip = r(:)' * y(:);
D = s * (1 - ip^2 / snr);

if(doDerivative)
    dr = K;
    dradj = Kadj;
    dD = 2 * s * Kadj(-ip * y / snr + r * (ip / snr)^2)';
    d2psi = 1;
end
end
