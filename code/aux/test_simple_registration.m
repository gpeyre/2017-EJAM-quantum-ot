%%
% Simple test for kernel-parameterized deformation.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_geometry/');
addpath('toolbox_quantum/tensor_logexp/');

sigma = .2; % width of the kernel
[K,Jac,push_fwd] = load_deformation_kernel(sigma);

% sampling points
p = 10; P = p*p;
z = linspace(-.3,1.3,p);
[y,x] = meshgrid(z,z);
Z = [x(:), y(:)];

% grid point for warping evaluation
n = 15; N = n*n;
x = linspace(0,1,n);
[y,x] = meshgrid(x,x);
X = [x(:), y(:)];

% test for a simple warping
c = [1/2 1/2];
R = sqrt( (X(:,1)-c(1)).^2 + ( X(:,2)-c(2) ).^2 );
V = .5 * [ -R.^2 .* (X(:,1)-c(1)), sqrt(R).* (X(:,2)-c(2))  ]; % displacement field 

% contraction 
phi = @(a,b)[-.22*sign(a).*abs(a).^2, ...
             -.22*sign(b).*abs(b).^2];
% contraction/expansion
phi = @(a,b)[-.22*sign(a).*abs(a).^2, ...
             +.50*sign(b).*abs(b).^2];         
% generic
phi = @(a,b)[-.05*b.^3, ...
             .4*a.^2];
                 
V =  phi(2*X(:,1)-1, 2*X(:,2)-1); % displacement field 

clf; hold on;
opt.edge_color = [0 0 0];
display_grid(  reshape(X + V, [n n 2]), opt );

% approximate using RBF
% |K*a-V|^2 + lambda*|a|^2
lambda = .1;
a = (K(X,Z)'*K(X,Z)+lambda*eye(P)) \ (K(X,Z)'*V);
V1 = K(X,Z)*a;

clf; hold on;
opt.edge_color = [0 0 0];
display_grid(  reshape(X + V, [n n 2]), opt );
opt.edge_color = [1 0 0];
display_grid(  reshape(X + V1, [n n 2]), opt );

aniso = .3;
e1 = cat(3,ones(n),zeros(n));
e2 = cat(3,zeros(n),ones(n));
mu = tensor_eigenrecomp(e1,e2,ones(n),aniso*ones(n));
mu = reshape(mu, [2 2 N]);

resh = @(mu)reshape(mu,[2 2 n n]);


% clf; hold on;
% plot_tensors_2d(  resh(mu), opt);
% display_grid(  reshape(X + V1, [n n 2]), opt );

%[e1,e2,l1,l2] = tensor_eigendecomp(resh(nu));
%nu1 = tensor_eigenrecomp(e1,e2,ones(n),.1*ones(n));

nu = push_fwd(mu,X,Z,a);

clf; hold on;
opt.scaling = 1.3/n;
plot_tensors_2d_scattered(reshape(nu, [2 2 N]), X+V1, opt);
opt.edge_color = [0 0 0];
display_grid(  reshape(X + V1, [n n 2]), opt );
 
% Tensor conjugacy   D^{-1} * T * D^{-1,T}






