%%
% Test for "transport-less" interpolation, i.e. just interpolation between
% two tensors. 

name = 'rot_small';
name = 'rot_large';
name = 'aniso-iso';
name = 'scaling';
name = 'ortho';
name = 'very_aniso';

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

barymode = 'soft';
barymode = 'hard';


rep = ['results/pointwise-interp/'];
[~,~] = mkdir(rep);

% generation of tensor from angle/aniso/scale
tensor = @(t,r,s)tensor_mult(tensor_creation(r, t), tensor_diag(s,s) );

aniso = .95;
switch name
    case 'rot_small'
        Q = tensor(0,aniso,1);
        R = tensor(.5*pi,aniso,1);
    case 'rot_large'
        Q = tensor(0,aniso,1);
        R = tensor(.9*pi,aniso,1);
    case 'very_aniso'
        aniso = 1-1e-9;
        Q = tensor(0,aniso,1);
        R = tensor(.9*pi,aniso,1);
    case 'ortho'
        Q = tensor(0,aniso,1);
        R = tensor(pi,aniso,1);
    case 'aniso-iso'
        Q = tensor(0,aniso,1);
        R = tensor(0,0,1);        
    case 'scaling'
        Q = tensor(0,0,1);
        R = tensor(0,0,.2);
end

tlist = linspace(0,1,15);
options.niter = 100;
options.barymode = barymode;
[P,err] = quantum_interp(Q,R, tlist, options);

opt.color = 'interp';
clf;
plot_tensors_1d(P, opt);
saveas(gcf, [rep 'interp-' name '-' barymode '.png'], 'png');
