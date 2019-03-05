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
function tests = SSD_opTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function oneDimTest(testCase)

% Create template image.
m = [99, 1];
Tc = ones(m);

% Create reference image.
Rc = ones(m);

% Set up domain.
omega = [0, 1];

% Create identity operator and its adjoint.
K = @(x) x;
Kadj = K;

% Compute SSD.
[Dc, r, dD, dr, dr_adj, d2psi] = SSD_op(Tc, Rc, omega, m, K, Kadj, 'doDerivative', true);
verifyEqual(testCase, Dc, 0);
verifyEqual(testCase, r, zeros(m));
verifyEqual(testCase, dD, zeros(m)');
verifyEqual(testCase, dr, K);
verifyEqual(testCase, dr_adj, Kadj);
verifyEqual(testCase, d2psi, prod(1 ./ m));

% Compute SSD for given scalar weight.
weights = 5;
[Dc, r, dD, dr, dr_adj, d2psi] = SSD_op(Tc, Rc, omega, m, K, Kadj, 'doDerivative', true, 'weights', weights);
verifyEqual(testCase, Dc, 0);
verifyEqual(testCase, r, zeros(m));
verifyEqual(testCase, dD, zeros(m)');
verifyEqual(testCase, dr, K);
verifyEqual(testCase, dr_adj, Kadj);
verifyEqual(testCase, d2psi, weights * prod(1 ./ m));

% Compute SSD for given matrix weight.
weights = 10 * speye(m(1), m(1));
[Dc, r, dD, dr, dr_adj, d2psi] = SSD_op(Tc, Rc, omega, m, K, Kadj, 'doDerivative', true, 'weights', weights);
verifyEqual(testCase, Dc, 0);
verifyEqual(testCase, r, zeros(m));
verifyEqual(testCase, dD, zeros(m)');
verifyEqual(testCase, dr, K);
verifyEqual(testCase, dr_adj, Kadj);
verifyEqual(testCase, d2psi, weights * prod(1 ./ m));

end

function oneDimDownsamplingTest(testCase)

% Create template image.
n = [10, 1];
Tc = ones(n);

% Create reference image.
m = [5, 1];
Rc = ones(m);

% Set up domain.
omega = [0, 1];

% Create downsampling operator and its adjoint.
K = @(x) 0.5 * (x(1:2:end) + x(2:2:end));
Kadj = @(y) 0.5 * kron(y, [1; 1]);

% Verify adjoint.
x = rand(n);
y = rand(m);
assert((K(x)'*y - x'*Kadj(y)) <= eps);

% Compute SSD.
[Dc, r, dD, dr, dr_adj, d2psi] = SSD_op(Tc, Rc, omega, m, K, Kadj, 'doDerivative', true);
verifyEqual(testCase, Dc, 0);
verifyEqual(testCase, r, zeros(m));
verifyEqual(testCase, dD, zeros(n)');
verifyEqual(testCase, dr, K);
verifyEqual(testCase, dr_adj, Kadj);
verifyEqual(testCase, d2psi, prod(1 ./ m));

end

function twoDimTest(testCase)

% Create template image.
m = [24, 31];
Tc = ones(m);

% Create reference image.
Rc = ones(m);

% Set up domain.
omega = [0, 1, 0, 1];

% Create identity operator and its adjoint.
K = @(x) x;
Kadj = K;

% Compute SSD.
[Dc, r, dD, dr, dr_adj, d2psi] = SSD_op(Tc, Rc, omega, m, K, Kadj, 'doDerivative', true);
verifyEqual(testCase, Dc, 0);
verifyEqual(testCase, r, zeros(m));
verifyEqual(testCase, dD, zeros(m)');
verifyEqual(testCase, dr, K);
verifyEqual(testCase, dr_adj, Kadj);
verifyEqual(testCase, d2psi, prod(1 ./ m));

end
