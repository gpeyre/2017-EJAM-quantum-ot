%%
% Test for Sinkhorn and barycenters on 2x2 matrices in a 2D domain.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('toolbox_anisotropic/');
addpath('data/images/');

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

name = '2d-smooth-rand';
name = '2d-iso-bump';
name = '2d-mixt-bump';
name = '2d-bump-donut';
name = '2d-aniso-fields';

rep = ['results/interpolation-2d/' name '/'];
[~,~] = mkdir(rep);

n = 32; % width of images
N = n*n; % #pixels
op = load_helpers(n);

opt.aniso = .06;
C = load_tensors_pair(name, n, opt);
mu = {}; 
for k=1:2    
    mu{k} = reshape(C{k},[2 2 N]);
end
n1 = 256; % upscaling for display
opt.laplacian = 'superbases';
opt.laplacian = 'fd';
opt.diffus_tau = .08;
opt.diffus_t = 50;
F = rendering_tensors_2d(mu,n1, [rep 'input'], opt);

switch name
    case '2d-aniso-fields'
        %% load a nice-looking rendering/coloring function
        render = @(x)texture_lut(x, 'red-metal');
    otherwise
        render = @(x)x;
end

%%
% Parameters

cost_type = 2;
if strcmp(name, '2d-aniso-fields')
    % cost_type = '2d-per';
end

% Ground cost
c = ground_cost(n,cost_type);
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
nu = quantum_interp(gamma, mu, m, 2, opt);
rendering_tensors_2d(nu,n1, [rep 'interpol']);

%%
% Compute an animation movie.

m = 60;
nu = quantum_interp(gamma, mu, m, 2, opt);
opt.disp_tensors = 0;
F = rendering_tensors_2d(nu,n1, '', opt);
Fr = {}; 
for k=1:m
    Fr{k} = render(F(:,:,k));
end
% display
k = 0; clf;
for k=1:m*5
    k = k+1;
    k1 = 1+mod(k-1,2*m-1);
    if k1>m
        k1=2*m-k1;
    end
    imageplot(Fr{k1}); drawnow;
end
% saveas video
opt.quality = 50;
write_video(Fr, [rep 'interpol'], 'mp4', opt);

