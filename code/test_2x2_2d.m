%%
% Test for Sinkhorn and barycenters on 2x2 matrices in a 2D domain.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('data/images/');

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

name = '2d-smooth-rand';
name = '2d-iso-bump';
name = '2d-mixt-bump';
name = '2d-bump-donut';

rep = ['results/interpolation-2d/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

n = 32; % width of images
N = n*n; % #pixels
op = load_helpers(n);

C = load_tensors_pair(name, n);
mu = {}; 
for k=1:2    
    mu{k} = reshape(C{k},[2 2 N]);
end
n1 = 256; % upscaling for display
rendering_tensors_2d(mu,n1, [rep 'input']);

%%
% Parameters

% Ground cost
c = ground_cost(n,2);
% regularization
epsilon = (.15)^2;  % large
epsilon = (.04)^2;  % small
epsilon = (.08)^2;  % medium
% fidelity
rho = 10;  %large
rho = 1;  %medium

%%
% Compute the coupling using Sinkhorn. 

options.niter = 500; % ok for .05^2
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 100;
opt.disp_tensors = 1;
nu = compute_quantum_interp(gamma, mu, m, 2, opt);
rendering_tensors_2d(nu,n1, [rep 'interpol']);

%%
% Compute an animation movie.


m = 60;
nu = compute_quantum_interp(gamma, mu, m, 2, opt);
opt.disp_tensors = 0;
F = rendering_tensors_2d(nu,n1, '', opt);
%
k = 0; clf;
while true
    k = k+1;
    k1 = 1+mod(k-1,2*m-1);
    if k1>m
        k1=2*m-k1;
    end
    imageplot(F(:,:,k1)); drawnow;
end

opt.quality = 50;
write_video(F, [rep 'interpol'], 'mp4', opt);

