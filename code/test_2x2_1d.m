%%
% Test for Sinkhorn and Barycenrters on 2x2 matrices in a 1D domain.

addpath('toolbox/');
addpath('toolbox_quantum/');

name = 'split';
name = 'multi-orient';
name = 'dirac-pairs';
name = 'iso-orient';
name = 'cross-orient';

trM = @(x)squeeze( x(1,1,:,:)+x(2,2,:,:) );

% mode of computation of logM/expM
global logexp_fast_mode;
logexp_fast_mode = 1;

rep = ['results/interpolation-1d/' name '/'];
if not(exist(rep))
    mkdir(rep); 
end    

% number of tensors in each measure
N = [1,1]*200;
N = [1,1]*60;
% number of barycenters
m = 8; 
% for display purposes
options.nb_ellipses = 30;

if strcmp(name, 'dirac-pairs')
N = [1,1]*25;
options.nb_ellipses = 25;
end

mu = load_tensors_pair(name,N);

% display
clf;
plot_tensors_1d(mu, options);
saveas(gcf,[rep 'input.eps'], 'epsc');
saveas(gcf,[rep 'input.png'], 'png');
% linear interplations
muL = {};
for k=1:m
	t = (k-1)/(m-1);    % TODO: correct inversion
	muL{k} = mu{1}*(1-t) + mu{2}*t;
end
clf;
plot_tensors_1d(muL, options);
saveas(gcf,[rep 'linear-interp.eps'], 'epsc');
saveas(gcf,[rep 'linear-interp.png'], 'png');

%%
% Parameters

% Ground cost
c = ground_cost(N,1);
% regularization
epsilon = (.15)^2;  % large
epsilon = (.08)^2;  % medium
epsilon = (.04)^2;  % small
% fidelity
rho = 10;  %large
rho = 1;  %medium
% prox param
lambda = rho/epsilon;



%%
% Run Sinkhorn.

options.niter = 30; % ok for .15^2
options.niter = 300; % ok for .05^2
options.disp_rate = 10;
options.tau = 1/(lambda+1); % should scale like 2/(1+lambda)

[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);
% trace of the couplings, to see where mass is flowing
G = squeeze( trM(gamma) );
G1 = G ./ repmat( sum(G,1), [size(G,1),1] );
G2 = G ./ repmat( sum(G,2), [1 size(G,2)] );
clf;
imageplot({G G1+G2 G1 G2},'',2,2);

%%
% Compute interpolation using an heuristic formula.

nu = quantum_interp(gamma, mu, m, 1);
% display 1D evolution as ellipses
clf;
plot_tensors_1d(nu, options);
saveas(gcf,[rep 'interp-ellipses.eps'], 'epsc');
saveas(gcf,[rep 'interp-ellipses.png'], 'png');

%%
% Barycenters.

options.tau = 1/(lambda+1); 
options.tau_nu = 1/options.tau; % not used anymore
% single test
options.niter = 120;
options.disp_rate = 5;
options.over_iterations = 1; % not needed anymore
w = [1/2 1/2];
options.tau = 1/(lambda+1); 
[nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);

clf; f = @(x)log10(x);
subplot(3,1,1); plot(f(err(:,1))); axis tight;
subplot(3,1,2); plot(f(err(:,2))); axis tight;
subplot(3,1,3); plot(f(err(:,3))); axis tight;


options.niter = 200; % sinkhorn #iterates
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
saveas(gcf,[rep 'barycenters-trace.eps'], 'epsc');
saveas(gcf,[rep 'barycenters-trace.png'], 'png');

% display 1D evolution as ellipses
clf;
plot_tensors_1d(nu, options);
saveas(gcf,[rep 'barycenters-ellipses.eps'], 'epsc');
saveas(gcf,[rep 'barycenters-ellipses.png'], 'png');