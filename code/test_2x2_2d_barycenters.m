%%
% Test for Sinkhorn and barycenters on 2x2 matrices in a 2D domain.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('data/images/');

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

name = '2d-bary';

rep = ['results/barycenters-2d/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

n = 12; % width of images
N = n*n; % #pixels
op = load_helpers(n);

C = load_tensors_pair(name, n);
mu = {}; 
for k=1:length(C)   
    mu{k} = reshape(C{k},[2 2 N]);
end

col = [1 0 0; 0 1 0; 0 0 1; 1 1 0];
opt.nb_ellipses = n;
for k=1:length(C)
    opt.color = col(k,:);
    clf;
    plot_tensors_2d(C{k}, opt);
    saveas(gcf, [rep 'input-' num2str(k) '.png'], 'png');
end

%%
% Ground cost

x = linspace(0,1,n);
[y,x] = meshgrid(x,x);
[X1,X2] = meshgrid(x(:),x(:));
[Y1,Y2] = meshgrid(y(:),y(:));
c0 = (X1-X2).^2 + (Y1-Y2).^2;
resh = @(x)reshape( x, [2 2 N N]);
flat = @(x)reshape(x, [2 2 N N]);
c = resh( tensor_diag(c0(:),c0(:)) );

%%
% Parameters

% regularization
epsilon = (.15)^2;  % large
epsilon = (.04)^2;  % small
epsilon = (.08)^2;  % medium
% fidelity
rho = 10;  %large
rho = 1;  %medium
% prox param
lambda = rho/epsilon;



%%
% Single test.

options.niter = 5;
options.disp_rate = NaN;
options.over_iterations = 10; % seems important to avoid oscilations
options.niter = 200; % sinkhorn #iterates
options.niter = 100;

m = 5; % number of barycenters
nu = {};
for k1=1:m
    for k2=1:m
        t1 = (k1-1)/(m-1);
        t2 = (k2-1)/(m-1);
        w = [(1-t1)*(1-t2), (1-t1)*t2 t1*(1-t2) t1*t2];
        fprintf('Barycenter (%d/%d,%d/%d):', k1, m, k2, m);
        [nu{k1,k2},gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
    end
end


for k1=1:m
    for k2=1:m
        t1 = (k1-1)/(m-1);
        t2 = (k2-1)/(m-1);
        w = [(1-t1)*(1-t2), (1-t1)*t2 t1*(1-t2) t1*t2];
        opt.color = sum( col .* repmat(w', [1 3]) );
        clf;
        plot_tensors_2d(reshape(nu{k1,k2}, [2 2 n n]), opt);
        drawnow;
        saveas(gcf, [rep 'barycenter-' num2str(k1) '-' num2str(k2) '.png'], 'png');
    end
end