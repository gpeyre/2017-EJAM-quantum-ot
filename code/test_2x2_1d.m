%%
% Test for Sinkhorn and Barycenrters on 2x2 matrices in a 1D domain.

addpath('toolbox/');
addpath('toolbox_quantum/');

name = 'split';
name = 'multi-orient';
name = 'iso-orient';
name = 'cross-orient';
name = 'dirac-pairs';
name = 'dirac-pairs-smooth';

trM = @(x)squeeze( x(1,1,:,:)+x(2,2,:,:) );

rep = ['results/interpolation-1d/' name '/'];
if not(exist(rep))
    mkdir(rep); 
end    

% number of tensors in each measure
N = [1,1]*200;
N = [1,1]*60;
N = [1,1]*30;
% number of barycenters
m = 9; 
% for display purposes
options.nb_ellipses = 30;
if strcmp(name, 'dirac-pairs')
    N = [1,1]*25;
    options.nb_ellipses = 25;
end

options.aniso = .95;
% options.aniso = 0; % un-comment this to force isotropic tensors.
mu = load_tensors_pair(name,N, options);

% display
clf;
plot_tensors_1d(mu, options);
saveas(gcf,[rep 'input.png'], 'png');
% linear interplations
muL = {};
for k=1:m
	t = (k-1)/(m-1);    % TODO: correct inversion
	muL{k} = mu{1}*(1-t) + mu{2}*t;
end
clf;
plot_tensors_1d(muL, options);
saveas(gcf,[rep 'linear-interp.png'], 'png');

%%
% Parameters

% mode of computation of logM/expM
global logexp_fast_mode;
logexp_fast_mode = 1;

% Ground cost
c = ground_cost(N,1);
% regularization
epsilon = (.15)^2;  % large
epsilon = (.04)^2;  % small
epsilon = (.01)^2;  % small
epsilon = (.06)^2;  % medium
% fidelity
rho = 1;  %medium
rho = 10;  %medium

%%
% Run Sinkhorn.

options.niter = 5000; 
options.disp_rate = 10;
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

% trace of the couplings, to see where mass is flowing
G = squeeze( trM(gamma) );
G1 = G ./ repmat( sum(G,1), [size(G,1),1] );
G2 = G ./ repmat( sum(G,2), [1 size(G,2)] );
clf;
imageplot({G G1+G2 G1 G2},'',2,2);

%%
% Compute interpolation using an heuristic formula.

options.sparse_mult = 30;
nu = quantum_interp(gamma, mu, m, 1, options);
% display 1D evolution as ellipses
clf;
plot_tensors_1d(nu, options);
saveas(gcf,[rep 'interp-ellipses.png'], 'png');

%%
% Test on a Barycenters.

options.tau = 1.8/(lambda+1); 
options.niter = 200; % sinkhorn #iterates
options.disp_rate = 5;
w = [1/2 1/2]; % weight, here iso-barycenter
options.tau = 1/(lambda+1); 
[nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);

% plot error decay
clf; f = @(x)log10(x);
subplot(3,1,1); plot(f(err(:,1))); axis tight;
subplot(3,1,2); plot(f(err(:,2))); axis tight;
subplot(3,1,3); plot(f(err(:,3))); axis tight;

%%
% All barycenters.

% full interpolation
nu = {};
for k=1:m
    t = (k-1)/(m-1);
    w = [1-t t];
    fprintf('Barycenter %d/%d:', k, m);
    [nu{k},gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
end

% display evolution of trace
clf; hold on;
for k=1:m
    t = (k-1)/(m-1);
    col = [t 0 1-t];
    plot(trM(nu{k}), 'color', col);
end
plot(trM(mu{1}), 'b--');
plot(trM(mu{2}), 'r--');
axis tight; box on;
saveas(gcf,[rep 'barycenters-trace.png'], 'png');

% display 1D evolution as ellipses
clf;
plot_tensors_1d(nu, options);
saveas(gcf,[rep 'barycenters-ellipses.png'], 'png');