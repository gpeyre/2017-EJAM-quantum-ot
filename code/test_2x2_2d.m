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

n = 15; % width of images
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
% prox param
lambda = rho/epsilon;
options.tau = 1/(lambda+1); 

%%
% Just compute the coupling using Sinkhorn. 

options.niter = 500; % ok for .05^2
options.disp_rate = NaN;
options.tau = 1/(lambda+1);
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

% Compute interpolation using an heuristic formula.
nu = compute_quantum_interp(gamma, mu, m, 2);
rendering_tensors_2d(nu,n1, [rep 'interpol']);


%%
% Single test.

options.niter = 5;
options.disp_rate = NaN;
options.over_iterations = 10; % seems important to avoid oscilations
w = [1/2 1/2];
% usual code, works fine


% [nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);

options.over_iterations = 1; % not needed anymore
options.niter = 50;
[nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options); 
rendering_tensors_2d(nu,n1, [], options);

f = @(x)log10(x);
clf;
subplot(3,1,1); plot(f(err(:,1))); axis tight;
subplot(3,1,2); plot(f(err(:,2))); axis tight;
subplot(3,1,3); plot(f(err(:,3))); axis tight;

%%
% full interpolation

options.niter = 200; % sinkhorn #iterates

m = 8; % numberof barycenters
nu = {};
for k=1:m
    t = (k-1)/(m-1);
    w = [1-t t];
    fprintf('Barycenter %d/%d:', k, m);
    [nu{k},gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
end
rendering_tensors_2d(nu,n1, [rep 'barycenters']);