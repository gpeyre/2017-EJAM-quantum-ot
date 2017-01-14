%% Compile

mex toolbox_quantum\tensor_logexp\tensorLog2x2.cpp ...
    -outdir toolbox_quantum\tensor_logexp -largeArrayDims ...
    -v COMPFLAGS="$COMPFLAGS /fp:fast /openmp" ...
    OPTIMFLAGS="$OPTIMFLAGS /Ox /openmp"
mex toolbox_quantum\tensor_logexp\tensorExp2x2.cpp ...
    -outdir toolbox_quantum\tensor_logexp -largeArrayDims ...
    -v COMPFLAGS="$COMPFLAGS /fp:fast /openmp" ...
    OPTIMFLAGS="$OPTIMFLAGS /Ox /openmp"

%%
% Test for Quantum-Sinkhorn, on a 1-D example.

clear

addpath('toolbox/');
addpath('toolbox_quantum/');

% name = 'split';
% name = 'multi-orient';
% name = 'cross-orient';
% name = 'iso-orient';
name = 'dirac-pairs';

% mode of computation of logM/expM
global logexp_fast_mode;
logexp_fast_mode = 4;

% number of tensors in each measure
N = [1,1]*200;
N = [1,1]*45;

% number of barycenters
m = 8; 

% for display purposes
options.nb_ellipses = 15;

if strcmp(name, 'dirac-pairs')
    N = [1,1]*25;
    options.nb_ellipses = 25;
end

mu = load_tensors_pair(name,N);

% display
clf;
plot_tensors_1d(mu, options);
title('Input (each row is a different distribution)');

% linear interplations
muL = {};
for k=1:m
	t = 1 - (k-1)/(m-1);    % TODO: correct inversion
	muL{k} = mu{1}*(1-t) + mu{2}*t;
end

figure;
plot_tensors_1d(muL, options);
title('Linear interpolation (y=t)');

%%
% Parameters

% regularization
% epsilon = (.15)^2;  % large
epsilon = (.08)^2;  % medium
% epsilon = (.04)^2;  % small

% fidelity
% rho = 10;  %large
rho = 1;  %medium

% prox param
lambda = rho/epsilon;

%%
% Ground cost

x = linspace(0,1,N(1));
y = linspace(0,1,N(2));
[Y,X] = meshgrid(y,x);
c0 = abs(X-Y).^2;
resh = @(x)reshape( x, [2 2 N(1) N(2)]);
flat = @(x)reshape(x, [2 2 N(1)*N(2)]);
c = resh( tensor_diag(c0(:),c0(:)) );

%%
% Barycenters.

options.niter = 300; % ok for .05^2
options.disp_rate = 10;
options.tau = .5 * 2/(lambda+1); % should scale like 2/(1+lambda)

options.tau = 1/(lambda+1); 
options.tau_nu = 1/options.tau; 

% single test
options.niter = 80;
options.disp_rate = 5;
options.over_iterations = 30; % seems important to avoid oscilations
w = [1/2 1/2];

% usual code, works fine
options.tau = 1/(lambda+1); 
options.disp_func = [];

global logexp_fast_mode;
logexp_fast_mode = 1;
tic
[nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
matlabTime = toc;

fprintf('Matlab time:  %g\n', matlabTime);

logexp_fast_mode = 4;
tic
[nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
cppTime = toc;

fprintf('C++ time:  %g\n', cppTime);

options.expMfunc = @tensorExp2x2;
optoins.logMfunc = @tensorLog2x2;

tic
[nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
defaultTime = toc;

fprintf('Default time:  %g\n', defaultTime);

%%

% new iterations, only works for small epsilon apparently :(
options.tau = 1/(lambda+1); 
options.tau_nu = .5/options.tau; 

options.niter = 800; % sinkhorn #iterates

% full interpolation
nu = {};
tic
for k=1:m
    t = (k-1)/(m-1);
    w = [t 1-t];
    fprintf('Barycenter %d/%d:', k, m);
    [nu{k},gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
end
toc

% display 1D evolution as ellipses
figure;
plot_tensors_1d(nu, options);
title('Transport');