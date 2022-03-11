function rec2=proxTV_toolbox(R,lambda)
%This script requires the toolbox_optim:
% matlab-toolboxes: Matlab toolboxes from www.numerical-tours.com.
% GitHub: https://github.com/gpeyre/matlab-toolboxes 

% This script computes prox_{lambda*TV}(R) via
% primal dual hybrid gradient method with Douchglas Rachford splitting
K=@(x) x;
Kadj=@(x) x;
K2 = @(x) grad(x);
KS2 = @(x) -div(x);
Amplitude = @(u) sqrt(sum(u.^2, 3));
F2 = @(u) lambda * sum(sum(Amplitude(u)));
Normalize = @(u) u ./ repmat(max(Amplitude(u), 1e-10), [1, 1, 2]);
ProxF2 = @(u, tau) repmat(perform_soft_thresholding(Amplitude(u), lambda * tau), [1 1 2]) .* Normalize(u);
ProxFS = compute_dual_prox(ProxF2);
F = @(x) norm(R - x, 'fro')^2/2;
G = @(x) 0;
ProxF = @(x,tau) (x - tau * R) / (1 + tau / 2);
ProxG = @(x,tau) x;
options.niter = 1000;
options.report = @(x) G(x) + F(K(x)) + F2(K2(x));


% Compute TV reconstruction with given template.
[rec2, ~] = perform_dr_pd(zeros(size(Kadj(R))), K, K2, Kadj, KS2, ProxF, ProxFS, ProxG, options);
end
function [x,R] = perform_dr_pd(x, K, K2, KS, KS2, ProxFS, ProxFS2, ProxG, options)

options.null = 0;
report = getoptions(options, 'report', @(x)0);
niter = getoptions(options, 'niter', 100);

if isnumeric(K)
    K = @(x)K*x;
end  
if isnumeric(KS)
    KS = @(x)KS*x;
end      

%%%% Douchglas Rachford parameters %%%%
sigma = getoptions(options, 'sigma', -1);
tau   = getoptions(options, 'tau', -1);
if sigma<0 || tau<0
    [L,~] = compute_operator_norm(@(x) KS(K(x)),randn(size(x)));
    [L2,~] = compute_operator_norm(@(x) KS2(K2(x)),randn(size(x)));
    sigma = .99/(L);
    sigma2 = .99/(L2);
    tau = 1;
end

y = K(x);
y2 = K2(x);
clear R;
lam_k = 1.5;

for i=1:niter
    R(i) = report(x);
    tmp_domain = KS(y);
    p1 = KS2(y2);
    tmp_domain = tmp_domain + p1;
    tmp_domain = x -(tau/2)*tmp_domain;
    p1 = ProxG(tmp_domain,tau);
    w1 = 2*p1 -x;
    tmp = y + (sigma/2) * K(w1);
    p2 = ProxFS(tmp, sigma);
    w2 = 2*p2 -y;
    tmp = y2 + (sigma2/2) * K2(w1);
    p22 = ProxFS2(tmp, sigma2);
    w22 = 2*p22 -y2;
    tmp_domain = KS(w2);
    z1= KS2(w22);
    tmp_domain = tmp_domain + z1;
    z1 = w1 - (tau/2)*tmp_domain;
    x = x + lam_k*(z1-p1);
    tmp_domain = 2*z1 - w1;
    z2 = w2 + (sigma/2)*K(tmp_domain);
    y = y + lam_k*(z2 -p2);
    z22 = w22 + (sigma2/2)*K2(tmp_domain);
    y2 = y2 + lam_k*(z22 -p22);
end
x = p1;
end