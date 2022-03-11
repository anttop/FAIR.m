%This script is a modification of the LDDMMobjFctn file of the Matlab-base toolbox
%LagLDDDM. Only minimal changes to the original were made, which are all
%marked by the flag $$CHANGE and be found in lines 109, 120, 134/135 and
%141/142.
% For details and license info concerning the original see
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM.
%==============================================================================

% function  [[Jc,para,dJ,H] = PALMobjFctn(T,Rc,omega,m,vRef,xc,omegaV,mV,N,vc,K,Kadj)
%
% Objective Function for LDDMM using Lagrangian PDE solver 
%
% computes J(vc) = D(T(yc(vc)),Rc), where
%
% vc       - is a velocity field (stationary or instationary)
% yc(vc)   = trafo(wc,xc) is obtained by tracing characteristics using RK4
%            scheme
% Tc       = T(yc) = imgModel(T,omega,yc),
% D(Tc,Rc) = distance(Tc,Rc,omega,m)

%Further, a regularizer term is given:
% S        = regularizer, e.g., S(uc) = 0.5*uc'*B*uc.
% 

% Input:
%   T      - data for template image, Tc = imgModel(T,omega,yc)
%   Rc     - reference image on grid, Rc = imgModel(R,omega,xc)
%   omega  - representation of computational domain
%   m      - discretization size
%   vRef   - reference velocity
%   xc     - discretization of Omega
%   omegaV - computational domain for velocity field (can be larger)
%   mV     - discretization size for velocities
%   N      - number of time steps for RK4 scheme
%   vc     -  current velocities
%   K      - a linear operator (optional)
%   Kadj   - adjoint of K (optional)
%
% Output:
%  Jc      - current function value J(vc)
%  para    - struct {Tc=T(y(vc)), Rc, omega, m, yc=y(vc,xc), Jc}, for plots
%  dJ      - gradient of  D(T(yc(vc)),Rc) with respect to vc
%  H       - second dervivate of the regularizer S with respect to vc
%==============================================================================
function [Jc,para,dJ,H] = PALMobjFctn(T,Rc,omega,m,vRef,xc,omegaV,mV,N,vc,K,Kadj)

para = struct([]);
if nargin == 0,
  help(mfilename);
  runMinimalExample;
  return;
elseif ~exist('vc','var') || isempty(vc),
  % if wc is not an input argument, reports status
  if nargout == 1, Jc = 'LDDMM';  return; end;
  % report current settings
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  vc      = trafo('w0');
  nt      = regularizer('get','nt');
  fprintf('Large Deformation Diffeomorphic Mapping (LDDMM):\n');
  fprintf('   J(vc) = D(T(yc(vc)),Rc) + S(vc-vRef) != min\n');
  fprintf('  %20s : %s\n','m',dimstr(m));
  fprintf('  %20s : %s\n','omega',dimstr(omega));
  fprintf('  %20s : %s\n','IMAGE MODEL',imgModel);
  fprintf('  %20s : %s\n','DISTANCE',distance);
  fprintf('  %20s : %s\n','TRAFO',trafo);
  fprintf('  %20s : %s\n','#timeSteps',num2str(N));
  fprintf('  %20s : %s\n','nt',num2str(nt));
  fprintf('  %20s : %s\n','mV',dimstr(mV));
  fprintf('  %20s : %s\n','omegaV',dimstr(omegaV));
  fprintf('  %20s : %s\n','length(vc)',num2str(length(vc)));
  Jc = vc; % return starting guess
  return;
end;
tspan = [1 0];
dim   = numel(omega)/2;
nt    = round(numel(vc)/(dim*prod(m)))-1;

% do the work ------------------------------------------------------------
matrixFree   = regularizer('get','matrixFree');
doDerivative = (nargout>2);            % flag for necessity of derivatives

% check if operator and adjoint are present
operator = exist('K', 'var') && exist('Kadj', 'var');

% compute transformation, distance, and regularization and combine these
if nt<1
    [yc,dy,pTrafo] = getTrafoFromVelocityRK4(vc,xc,'omega',omegaV,'m',mV,'N',N,'tspan',tspan,'doDerivative',doDerivative);
else
    [yc,dy,pTrafo] = getTrafoFromInstationaryVelocityRK4(vc,xc,'omega',omegaV,...
                        'm',mV,'N',N,'tspan',tspan,'doDerivative',doDerivative);
end
[Tc,dT] = imgModel(T,omega,center(yc,m),'doDerivative',doDerivative);

% compute distance
if(operator)
    dist = str2func(distance());
    [Dc,~,dD,dres,dres_adj,d2psi] = dist(Tc,Rc,omega,m,K,Kadj,'doDerivative',true);
else
    [Dc,~,dD,dres,d2psi] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);
end

% compute regularizer
[Sc,dS,d2S] = regularizer(vc-vRef,omegaV,mV,'doDerivative',doDerivative,'tspan',tspan);

Jc = Dc; %+ Sc;   %CHANGE

% collect variables for plots
if(operator)
    Dshow = @(x,y,omega,m) viewImage(abs(K(x(:))-y));
    para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',yc,'Jc',Jc,'Sc',Sc,'Dc',Dc,...
        'N',pTrafo.N,'Dshow',Dshow);
else
    para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',yc,'Jc',Jc,'Sc',Sc,'Dc',Dc,...
        'N',pTrafo.N);
end
    
if ~doDerivative, return; end;
dD = dD*dT*dy;
%dJ = dD + dS;                   $$$ CHANGE   
dJ= dD;

if nargout<4, return; end;

% multiply outer and inner derivatives, note: dy might be sparse
if not(matrixFree),
   % dres = dres*dT*dy;             $$$ CHANGE
   % H = dres'*d2psi*dres + d2S;
     H= d2S;
else
    % include operator if present
    if(operator)
        dres    = @(x) dres(dT*dy*x);
        dres_adj= @(x) dy'*dT'*dres_adj(x);
    else
        dres_tmp = dres;
        dres = @(x) dres*dT*dy*x;
        dres_adj = @(x) dy'*dT'*dres_tmp'*x;
    end
    
    % approximation to d2D in matrix free mode
    % d2D   = dr'*d2psi*dr
    % P and P' are operators matrix free
    H.omega     = omegaV;
    H.m         = mV;
    H.nt        = nt;
    H.tspan     = tspan;
    H.d2D.how   = '*dr''*d2psi*dr';
    H.d2D.P     = @(x) x;
    H.d2D.dr    = dres;
    H.d2D.dr_adj= dres_adj;
    H.d2D.d2psi = d2psi;

    H.d2S = d2S;

end;

function runMinimalExample
setup2DGaussianData;
lvl = 5;
regularizer('reset','regularizer','mbElastic','alpha',1)
v0 = getVelocityStartingGuess(omega,m);
xc = getCellCenteredGrid(omega,ML{lvl}.m);
T  = reshape(imgModel(ML{lvl}.T,omega,center(xc,ML{lvl}.m)),ML{lvl}.m);
Rc = imgModel(ML{lvl}.R,omega,center(xc,ML{lvl}.m));

fctn = @(vc) LDDMMobjFctn(T,Rc,omega,ML{lvl}.m,v0,xc,omega,m,10,vc);
checkDerivative(fctn,v0)


