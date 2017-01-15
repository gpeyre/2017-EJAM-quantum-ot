%%
% Test for Sinkhorn convergence speed.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('toolbox_anisotropic/');
addpath('data/images/');

global logexp_fast_mode;
logexp_fast_mode = 4; % fast mex

name = '2d-bump-donut';

rep = 'results/speed/';
[~,~] = mkdir(rep);

n = 20; % width of images
N = n*n; % #pixels
op = load_helpers(n);

opt.aniso = .06;
mu = load_tensors_pair(name, n, opt);

resh = @(x)reshape(x, [size(x,1) size(x,2) size(x,3)*size(x,4)]);
iresh = @(x)reshape(x, [size(x,1) size(x,2) sqrt(size(x,3)) sqrt(size(x,3))]);

%%
% Sinkhorn

% Ground cost
c = ground_cost(n,2);
% regularization
epsilon = (.08)^2;  % medium
% fidelity
rho = 1;  %medium

options.niter = 500; 
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
[gamma,u,v,err] = quantum_sinkhorn(resh(mu{1}),resh(mu{2}),c,epsilon,rho, options);

%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 100;
opt.disp_tensors = 1;
fprintf('Interpolating: ');
nu = quantum_interp(gamma, {resh(mu{1}) resh(mu{2})}, m, 2, opt);

opt.nb_ellipses = n/2;
opt.scaling = .7;
for k=1:m
    t = (k-1)/(m-1);
    opt.color = [t 0 1-t];
    clf; plot_tensors_2d(iresh(nu{k}), opt); drawnow;
end

%% 
% Compute convergence rate of Sinkhorn for different relaxation parameter
% rho. 

options.niter = 250; 
alist = linspace(.8,1.99,20);
Err = [];
for k=1:length(alist)
    options.tau = alist(k)*epsilon/(rho+epsilon); 
    [gamma,u,v,err] = quantum_sinkhorn(resh(mu{1}),resh(mu{2}),c,epsilon,rho, options);
    Err(1:size(err,1),k) = err(:,1);
end

%%
% Display convergence.

clf;  hold on;
for k=1:length(alist)
    t = (k-1)/(length(alist)-1);
    plot(log10(Err(:,k)/max(Err(1,:))), 'color', [t 0 1-t]);
end
set(gca, 'FontSize', 15); 
axis([1 options.niter min(log10(Err(:))) 0]);
box on;
saveas(gcf, [rep 'convergence-curve.eps'], 'epsc');

%%
% Display Linear convergence rate.

i = options.niter; j = 50;
Q = - ( log(Err(i,:))-log(Err(j,:)) ) / (i-j);
clf; hold on;
for k=1:length(alist)-1
    t = (k-1)/(length(alist)-2);
    plot(alist(k:k+1), Q(k:k+1), 'LineWidth', 2, 'color', [t 0 1-t]); 
end
set(gca, 'FontSize', 15); 
axis([min(alist) max(alist) min(Q) max(Q)]); box on;
saveas(gcf, [rep 'convergence-rate.eps'], 'epsc');



