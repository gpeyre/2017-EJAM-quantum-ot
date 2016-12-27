%%
% Test for Sinkhorn and barycenters on 2x2 matrices in a 2D domain, using sparse
% coupling computation to handle large scale problems.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('data/images/');

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

name = '2d-bump-donut';
name = 'struct-tensors';

rep = ['results/interpolation-2d-sparse/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

nsub = 3; % number of coarsening steps
p = 16; % low res-image width
P = p*p;
n0 = p*2^nsub; % high res-image width
N0 = n0*n0; % #pixels

switch name
    case 'struct-tensors'
        name_img = {'hibiscus', 'flower'};
        op = load_helpers(n0);
        sigma = 1.5; % width of averaging for the ST
        for k=1:2
            f0{k} = load_image(name_img{k}, n0); 
            f0{k} = rescale(sum(f0{k},3));
            T0{k} = op.ST(f0{k},sigma);
            C{k} = op.T2C(T0{k});
            mu0{k} = reshape(C{k},[2 2 N0]);
            % display
            clf;  opt.sub = 4;
            plot_tensor_field(op.rotate(T0{k}), f0{k}, opt);
            saveas(gcf, [rep 'input-images-' num2str(k) '.png'], 'png')
        end
    otherwise
        C = load_tensors_pair(name, n0);
        mu0 = {};
        for k=1:2
            mu0{k} = reshape(C{k},[2 2 N0]);
        end
end
n_render = 256; % upscaling for display

mu = tensor_coarsening(mu0,2,nsub);
% rendering_tensors_2d(mu,n_render, [rep 'input']);

%%
% Ground cost

x = (0:1/p:1-1/p);
[y,x] = meshgrid(x,x);
[X1,X2] = meshgrid(x(:),x(:));
[Y1,Y2] = meshgrid(y(:),y(:));
c0 = (X1-X2).^2 + (Y1-Y2).^2;
resh = @(x)reshape(x, [2 2 P P]);
c = resh( tensor_diag(c0(:),c0(:)) );

%%
% Parameters

% regularization
epsilon = (.15)^2;  % large
epsilon = (.08)^2;  % medium
epsilon = (.06)^2;  % small
% fidelity
rho = 10;  %large
rho = 1;  %medium


%%
% Full sinkhorn

options.niter = 2*4000; % ok for .05^2
options.disp_rate = NaN;
options.tau = 1.8 * epsilon/(rho+epsilon); % use extrapolation to speed-up
options.tol = 1e-13;
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

% Compute McCann interpolation.
m = 9; % number of barycenters
opt.sparse_mult = 10; % sparsification of the coupling, to speed up
nu = compute_quantum_interp(gamma, mu, m, 2, opt);
rendering_tensors_2d(nu{(m+1)/2},n_render, [rep 'isobary']);

% sparse cost
E = trM(gamma,1);
if 1
    sparsity = 6; % sparsity factor
    v = sort(E(:), 'descend');  v = v(round(P*sparsity)); % threshold
else
    v = 1e-5;
end
cS = sparse(c0 .* (E>v));

% same, but on a sparse grid
[gammaS,u,v,err2] = quantum_sinkhorn(mu{1},mu{2},cS,epsilon,rho, options);
nu = compute_quantum_interp(gammaS, mu, m, 2);
rendering_tensors_2d(nu{(m+1)/2},n_render, [rep 'isobary-sparse']);

%%
% Multi-scale sparse sinkhorn.
 
[i,j,~] = find(cS); % current grid
for i_sub=1:nsub  
    % image width at this scale
    n = p*2^i_sub;
    % up-sample by factor 2 the grid.   
    [i,j,cS] = grid_subdivide(i,j,n/2, 2, 1);   
    mu1 = tensor_coarsening(mu0,2,nsub-i_sub);
    % sparse sinkhorn
    [gammaS,u,v,err2] = quantum_sinkhorn(mu1{1},mu1{2},cS,epsilon,rho, options);
    % sparsify grid
    E = trM(gammaS.T,1);
    v = sort(E(:), 'descend');  v = v( round(n*n*sparsity) ); % threshold   
    i = i(E>v); j = j(E>v); 
end

nu = compute_quantum_interp(gammaS, mu1, m, 2);
rendering_tensors_2d(nu{(m+1)/2},n_render, [rep 'isobary-sparse']);

for k=1:m
    A = reshape( trM(nu{k}, 1), [n n] );
    imwrite(rescale(A), [rep 'barycenter-trace-' num2str(k) '.png'], 'png');
end

if strcmp(name, 'struct-tensors')
    % least square tensor reconstruction     
    for k=1:m
        t = (k-1)/(m-1);
        T = op.C2T( reshape(nu{k},[2 2 n n]) ); % target tensor
        % initialization : use linear interpolation
        X = (1-t)*sort(f0{1}(:)) + t * sort(f0{2}(:)); % target historagm using OT
        f_init = (1-t)*f0{1} + t * f0{2};
        f_init = perform_histogram_equalization(f_init,X);
        opt.lambda = 2.5*1e-5;
        opt.niter_bfgs = 200;
        f1 = struct_tensor_reconstruc(f_init, T, sigma, opt);
        f1 = perform_histogram_equalization(f1,X);
        fR{k} = f1;
        imwrite(rescale(f1), [rep 'barycenter-reconstr-' num2str(k) '.png'], 'png');
    end
end
