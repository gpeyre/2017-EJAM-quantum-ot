%%
% Test for Sinkhorn on 2x2 matrices in a 2D domain.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('toolbox_anisotropic/');
addpath('data/images/');

name = '2d-smooth-rand';
name = '2d-iso-bump';
name = '2d-mixt-bump';
name = '2d-aniso-fields';
name = '2d-bump-donut';

rep = ['results/interpolation-2d/' name '/'];
[~,~] = mkdir(rep);

n = 50; % width of images
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

if 0
    % just to display the colorbar
    A = repmat( linspace(0,1,256)', [1 25] ); 
    A = opt.render(A);
    imwrite(A, [rep 'colorbar.png'], 'png');
end

%%
% Compute the coupling using Sinkhorn. 


global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

% Ground cost
c = ground_cost(n,2);
% regularization
epsilon = (.08)^2;  % medium
% fidelity
rho = 1;  %medium
% run sinkhorn
options.niter = 500; % ok for .05^2
options.disp_rate = NaN; % no display
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
fprintf('Sinkhorn: ');
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 100;
opt.disp_tensors = 1;
fprintf('Interpolating: ');
nu = quantum_interp(gamma, mu, m, 2, opt);

%%
% Display with a background texture.

[F,Fr] = rendering_tensors_2d(nu,n1, [rep 'interpol'], opt);
for k=1:m
    imwrite( Fr{k}, [rep 'interpol-render-' num2str(k) '.png'], 'png' );
end

%%
% Display with no backgound texture.

opt.nb_ellipses = n/2;
opt.scaling = .7;
for k=1:m
    t = (k-1)/(m-1);
    opt.color = [t 0 1-t];
    clf; plot_tensors_2d(reshape(nu{k}, [2 2 n n]), opt); drawnow;
    saveas(gcf, [rep 'interpol-ellipses-' num2str(k) '.png']);
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

