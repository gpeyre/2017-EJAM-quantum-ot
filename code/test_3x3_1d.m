%%
% Test for Sinkhorn and barycenters on 3x3 matrices in a 1D domain.


addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('data/images/');

name = 'plate-elong';

rep = ['results/interpolation-3d/' name '/'];
if not(exist(rep))
    mkdir(rep); 
end    


global logexp_fast_mode;
logexp_fast_mode = 1;

% number of tensors
n = 60;
% number of barycenters
m = 8;

%%
% Create tensor fields

mu = load_tensors_pair(name,[n n]);

opt.nb_ellipses = min(n,30);
clf;
plot_tensors_volume(mu, opt);

% linear interpolation
muL = {};
for k=1:m
	t = (k-1)/(m-1);
	muL{k} = mu{1}*(1-t) + mu{2}*t;
end
clf;
plot_tensors_volume(muL, opt);
saveas(gcf,[rep 'linear-interp.png'], 'png');


%%
% Parameters

% regularization
epsilon = (.04)^2;  % small
epsilon = (.15)^2;  % large
epsilon = (.06)^2;  % medium
% fidelity
rho = 10;  %large
rho = 1;  %medium
% prox param
lambda = rho/epsilon;

%%
% Ground cost

x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);
c0 = abs(X-Y).^2;
resh = @(x)reshape( x, [1 1 n n] );
c = zeros(3,3,n,n);
for k=1:3
    c(k,k,:,:) = resh(c0);
end

%%
% Run Sinkhorn.

options.niter = 600; % ok for .05^2
options.disp_rate = NaN;
options.tau = 1/(lambda+1); % should scale like 2/(1+lambda)
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);


%%
% Compute interpolation using an heuristic formula.

nu = quantum_interp(gamma, mu, m, 1);
% display 1D evolution as ellipses
clf;
plot_tensors_volume(nu, opt);
saveas(gcf,[rep 'interp-ellipses.png'], 'png');




%%
% Barycenters.

options.tau = 1/(lambda+1); 
% single test
options.niter = 250;

% full interpolation
nu = {};
for k=1:m
    t = (k-1)/(m-1);
    w = [t 1-t];
    fprintf('Barycenter %d/%d:', k, m);
    [nu{k},gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
end


% display 1D evolution as ellipses
clf;
plot_tensors_volume(nu, opt);
saveas(gcf,[rep 'barycenters-ellipses.png'], 'png');