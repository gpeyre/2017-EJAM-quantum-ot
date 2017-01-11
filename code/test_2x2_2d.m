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
name = '2d-aniso-fields';
name = '2d-bump-donut';

rep = ['results/interpolation-2d/' name '/'];
[~,~] = mkdir(rep);

n = 32; % width of images
N = n*n; % #pixels
op = load_helpers(n);

opt.aniso = .06;
C = load_tensors_pair(name, n, opt);


opt.diffus_tau = .08;
opt.diffus_t = 50;
if strcmp(name, '2d-bump-donut')
    for k=1:2
        [e1,e2,l1,l2] = tensor_eigendecomp(C{k});
        l1 = max(l1,.03);
        l2 = max(l2,.03);
        C{k} = tensor_eigenrecomp(e1,e2,l1,l2);        
    end
    opt.diffus_t = 20;
end

mu = {}; 
for k=1:2    
    mu{k} = reshape(C{k},[2 2 N]);
end
n1 = 256; % upscaling for display
opt.laplacian = 'superbases';
opt.laplacian = 'fd';
opt.disp_tensors = 1;
opt.render = @(x)texture_lut(x, 'red-metal');
[F,Fr] = rendering_tensors_2d(mu,n1, [rep 'input'], opt);

%%
% Parameters

cost_type = 2;
if strcmp(name, '2d-aniso-fields')
    % cost_type = '2d-per'; %TODO: FIX THIS
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
[F,Fr] = rendering_tensors_2d(nu,n1, [rep 'interpol'], opt);

for k=1:m
    imwrite( Fr{k}, [rep 'interpol-render-' num2str(k) '.png'], 'png' );
end
%%
% Compute an animation movie.

m = 60;
nu = quantum_interp(gamma, mu, m, 2, opt);

if 0 % strcmp(name, '2d-aniso-fields')
    % remap anisotropy
    for k=1:m
        [e1,e2,~,~] = tensor_eigendecomp(reshape(nu{k}, [2 2 n n]));
        nu{k} = tensor_eigenrecomp(e1,e2,ones(n),ones(n)*opt.aniso);
        nu{k} = reshape(nu{k}, [2 2 n*n]);
    end
end

opt.disp_tensors = 0;
[F,Fr] = rendering_tensors_2d(nu,n1, '', opt);
% display
k = 0; clf;
for k=1:m
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

