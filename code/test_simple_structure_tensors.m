%%
% Test for structure tensors.

name = 'hibiscus';

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('data/images/');

%%
% Image loading.

n = 256;
f0 = load_image(name, n); 
f0 = rescale(sum(f0,3));

%%
% Helpers. 

op = load_helpers(n);

% display an example of structure tensor. 
options.sub = 8;
clf; sigma0 = 2;
plot_tensor_field(op.rotate(op.ST(f0,sigma0)), f0, options);
title(['\sigma=' num2str(sigma0)]);

%%
% Fitting energy and its gradient. 


% check gradient validity
sigma0 = 1;
T0 = op.ST(f0,sigma0);
func = @(f)deal( op.TensorFit.E(f,T0,sigma0), op.TensorFit.nablaE(f,T0,sigma0) );
[e,deriv,deriv_fd] = check_gradient(func, [n n], 10);
clf; hold on;
plot(deriv, 'b.', 'MarkerSize', 10);
plot(deriv_fd, 'ro');

%% 
% L-BFGS parameters

flat = @(x)x(:); resh = @(x)reshape(x, [n n]);
BfgsCallback = @(x,t,sigma)deal( op.TensorFit.E(resh(x),t,sigma), ...
    flat(op.TensorFit.nablaE(resh(x),t,sigma)) );
opt.niter = 100;
opt.bfgs_memory = 5; 
opt.report = @(x,v)v;

%%
% Denoising by comparing to the original tensors

sigma0 = 1;
T0 = op.ST(f0,sigma0);
fInit = f0 + randn(n,n)*.06;
[f, R, info] = perform_bfgs(@(x)BfgsCallback(x,T0,sigma0), fInit(:), opt);
f = resh(f);
clf; plot(R); 
clf; imageplot({fInit f});

%%
% Sharpening by increasing ansotropy

[e1,e2,l1,l2] = tensor_eigendecomp(op.T2C(T0));
[a,e] = eigen_remaper(l1,l2,+1);
[L1,L2] = eigen_remaper(a.^(.4),e,-1);
C1 = tensor_eigenrecomp(e1,e2,L1,L2);
T1 = op.C2T(C1);

clf; sigma0 = 2;
plot_tensor_field(op.rotate(T1), f0, options);

[f, R, info] = perform_bfgs(@(x)BfgsCallback(x,T1,sigma0), f0(:), opt);
f = resh(f);

clf; plot(R); 
clf; imageplot({f0 f});
        