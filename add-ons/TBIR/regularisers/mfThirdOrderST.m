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
%
% function [Sc, dS, d2S] = mfThirdOrderST(uc,omega,m,varargin)
%
% Matrix-free spatio-temporal third-order regularization energy for vc
% where vc is cell-centered
%
% Sv) = 0.5 * \int_{\omega}\int_0^1 
%               alpha(1)*v(x,t)'*A*v(x,t) + alpha(2)*v(x,t)'*B*v(x,t)
%               + alpha(3)*v(x,t)'*v(x,t) dx dt,
%
% where A is the third-order partial derivative operator and B the
% first-order time derivative operator.
%
% Input:
%
%   vc          instationary velocity field (cell-centered)
%   omega       spatial domain
%   m           number of discretization points in space
%   varargin    optional parameters (see below)
%
% Optional Input:
%   tspan       time span (default: [0 1])
%   nt          number of time points for velocity
%
%
% Output:
%
%   Sc          current value  (0.5 * hd * uc'* A *uc)
%   dS          derivative     (hd * uc'*A )
%   d2S         Hessian        A
%
%  if ~matrixFree,  d2S is sparse matrix; else, d2S is struct; endif
function [Sc, dS, d2S] = mfThirdOrderST(vc, omega, m, varargin)
    if nargin == 0
        help(mfilename);
        return;
    end

    if strcmp(vc,'para')
        Sc = 'cell-centered';       % grid
        dS = 1;                     % matrixFree
        d2S = @spectralPrecondPCG;  % solver
        return;
    end

    if ~exist('mOld','var'),     mOld = [];     end
    if ~exist('omegaOld','var'), omegaOld = []; end
    if ~exist('alphaOld','var'), alphaOld = []; end

    alpha       = [1 1e-3];
    tspan       = [0 1];
    nt          = [];
    for k=1:2:length(varargin) % overwrites default parameter
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end

    dim = numel(omega)/2;

    if isempty(nt) % roughly estimate nt
        nt = round(numel(vc)/(prod(m)*dim))-1;
    end

    d2S.regularizer = regularizer;
    d2S.alpha  = alpha;
    d2S.B      = @(omega,m)     getDerivativeMatrixST(omega,tspan,m,nt,alpha);
    d2S.d2S    = @(u,omega,m)   derivativeOperatorST(u,omega,tspan,m,nt,alpha);
    d2S.diag   = @(omega,m)     getDerivartiveDiag(omega,tspan,m,nt,alpha); %getDiag(omega,tspan,m,nt,alpha);
    d2S.solver = @spectralPrecondPCG;
    d2S.res    = vc;
    dS         = d2S.d2S(vc,omega,m)';
    Sc         = 0.5*dS*vc;
end

function getDerivativeMatrixST(omega,tspan,m,nt,alpha)
    error('Not implemented!');
end

% matrix free implementation of spatio-temporal derivative operator
function Ay = derivativeOperatorST(uc,omega,tspan,m,nt,alpha)
    dim = numel(omega)/2;
    h   = (omega(2:2:end)-omega(1:2:end))./m;
    hd  = prod(h(1:dim));
    dt  = abs(tspan(1)-tspan(2))/nt;
    w   = dt*[1/2;ones(nt-1,1);1/2];
    n   = prod(m);

    uc   = reshape(uc,dim*n,[]);
    ntp1 = size(uc,2);
    Ay = 0*uc;

    switch dim
        case 2
            D2x =  D2(1,omega,m);
            D2y =  D2(2,omega,m);

            %%% Compute third-order derivatives.
            for k=1:ntp1
                % Define first derivatives and transposes.
                d1 = @(Y) (Y(2:end,:)-Y(1:end-1,:))/h(1);
                d2 = @(Y) (Y(:,2:end)-Y(:,1:end-1))/h(2);

                d1T = @(Y) ([-Y(1,:);Y(1:end-1,:)-Y(2:end,:);Y(end,:)])/h(1);
                d2T = @(Y) ([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)])/h(2);

                uct = reshape(uc(:,k),[m dim]);
                t1 = D2x*d1T(d1(D2x*uct(:, :, 1))) + 4*d2T(D2x*D2x*d2(uct(:, :, 1))) + 4*d1T(d1(uct(:, :, 1)*D2y))*D2y + d2T(d2(uct(:, :, 1)*D2y))*D2y;
                t2 = D2x*d1T(d1(D2x*uct(:, :, 2))) + 4*d2T(D2x*D2x*d2(uct(:, :, 2))) + 4*d1T(d1(uct(:, :, 2)*D2y))*D2y + d2T(d2(uct(:, :, 2)*D2y))*D2y;
                Ay(:,k) = (w(k)*hd*alpha(1))*[t1(:);t2(:)];
            end
            Ay = Ay(:);
        otherwise
            error('%s - dimension %d not supported.',mfilename,dim);
    end
    % time: diffusion
    if alpha(2) > 0
        d3  = @(Y) (Y(:,2:end)-Y(:,1:end-1))/dt;
        d3T = @(Y) ([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)])/dt;
        At = d3T(d3(uc));
        Ay = Ay(:) + (alpha(2)* hd*dt)*At(:);
    end
    % L2 regularization.
    if alpha(3) > 0
        Ay = Ay(:) + alpha(3) * hd*dt * uc(:);
    end
end

function D = getDerivartiveDiag(omega, tspan, m, nt, alpha)
    error('Not implemented!');
end

function D = D2(i,omega,m)
    h = (omega(2:2:end)-omega(1:2:end))./m;
    D = spdiags(ones(m(i),1)*[1, -2, 1], -1:1, m(i), m(i)) / h(i)^2;
    D([1,end]) = -D([2,end-1]);
end
