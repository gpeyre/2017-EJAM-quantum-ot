%%
% Test for anisotropic diffusion.

if exist('batch_mode_sigma') && batch_mode_sigma==1
    batch_mode=0;
    for sigma = 2:2:50
        fprintf('-------- sigma=%d ---------\n', sigma);
        test_diffusion;
    end    
end

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('data/images/');

rep = 'results/diffusion/';
if not(exist(rep))
    mkdir(rep);
end

n = 256;


%%
% Helpers. 

op = load_helpers(n);

%%
% Create random tensors.

name = '2d-smooth-rand';
name = '2d-iso-bump';
name = '2d-bump-donut';

if not(exist('sigma'))
    sigma = 50;
end
options.sigma = sigma; 
mu = load_tensors_pair(name,n, options);
T1 = op.C2T(mu{1});

% display using diffusion
tau = .06; % diffusion step size
t = 10; % final time
f = anisotropic_diffusion(T1, randn(n), t, tau);
clf; imageplot(f);
% imwrite(rescale(f), [rep name 'sigma-' num2str(sigma) '.png'], 'png');

% display as ellipsoids
opt.nb_ellipses = 16;
opt.color = [1 0 0];
opt.color_edge = [.7 0 0];
opt.scaling = .5;
opt.image = f;
clf;
plot_tensors_2d(op.T2C(T1), opt);
saveas(gcf, [rep name 'sigma-' num2str(sigma) '.png'], 'png');
