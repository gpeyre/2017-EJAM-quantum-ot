%%
% Test for registration with Wasserstein fidelity.

addpath('data/images/');
addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_geometry/');
addpath('toolbox_quantum/tensor_logexp/');

% sampling points
p = 20; P = p*p;
z = linspace(-.3,1.3,p);
[y,x] = meshgrid(z,z);
Z = [x(:), y(:)];

sigma = .05; % width of the kernel
sigma = .05; % width of the kernel
[K,Jac,push_fwd] = load_deformation_kernel(sigma);

% load input measures to register together. 

name = 'simple';
name = 'structure-tensors';

names = {'bar-vert-1', 'bar-vert-2'};
names = {'curve-1', 'curve-2'};


switch name
    case 'simple'
        % just for debuggin, simple case of pairs of diracs
        q = 10; % to dimensionate plots
        mu = {tensor_id([1;1], 2) tensor_id([1;1], 2)};
        a = .3; % spacing along y
        b = .8; % spacing along x
        X = { .5+[-b/2,-a/2; +b/2,-a/2], .5+[-b/2,+a/2; +b/2,+a/2] }; 
        N = [2 2];
    case 'structure-tensors'
        % structure tensors
        q = 40; % input size
        sigma0 = 2; % ST width
        op = load_helpers(q);
        % number of point in the measure
        N = round(3*[q q]);
        % input grid
        x0 = (0:q-1)'/q;
        [Y0,X0] = meshgrid(x0,x0);
        mu = {}; S = {}; f = {}; X = {};
        for k=1:2
            f{k} = load_image(names{k}, q);
            f{k} = rescale(sum(f{k},3));
            S{k} = op.rotate( op.ST(f{k},sigma0) );
            % detect large amplitude tensors
            a = S{k}(:,:,1)+S{k}(:,:,2);
            [~,I] = sort(a(:), 'descend'); I = I(1:N(k));
            X{k} = [X0(I) Y0(I)];
            mu{k} = zeros(2,2,N(k));
            v = S{k}(:,:,1); mu{k}(1,1,:) = reshape(v(I), [1 1 N(k)]);
            v = S{k}(:,:,2); mu{k}(2,2,:) = reshape(v(I), [1 1 N(k)]);
            v = S{k}(:,:,3); mu{k}(1,2,:) = reshape(v(I), [1 1 N(k)]);
            v = S{k}(:,:,3); mu{k}(2,1,:) = reshape(v(I), [1 1 N(k)]);
        end
end

%%
% Display input tensor fields.

clf; hold on;
opt.scaling = 1/q;
opt.color = [1 0 0];
plot_tensors_2d_scattered(mu{1}, X{1}, opt);
opt.color = [0 0 1];
plot_tensors_2d_scattered(mu{2}, X{2}, opt);

%%
% Sinkhorn parameters.

global logexp_fast_mode;
logexp_fast_mode = 4; % fast mex

% regularization
epsilon = (.08)^2;  % medium
epsilon = (.06)^2;  % medium
% fidelity
rho = 1;  %medium
rho = 3;  %medium
options.niter_sinkhorn = 3000; 

% initialize the gradient descent on a.
a = zeros(P,2);
tau = .01;
niter = 30;
Kmat = K(X{1},Z); % fixed
% underlying grid, just for display
n0 = 18;
x0 = linspace(0,1,n0); [Y0,X0] = meshgrid(x0,x0); X0 = [X0(:) Y0(:)];
use_a = 1; 
% for "simple"
tau = .3;
sigma_grad = 3; % introduce an extra regularization of the gradient on a
% use_a = 0; tau = .3; 
%
Y = X{1}; 
for k=1:niter
    %%% displace the input measure %%%
    if use_a
        Y = X{1} + Kmat*a;
        nu = push_fwd(mu{1}, X{1}, Z, a);
    else
        nu = mu{1};
    end
    Nu{k} = nu;
    [W,nabla_x,gamma] = quantum_fidelity(nu, mu{2}, Y', X{2}', rho, epsilon, options);
    %[W,nabla_x,gamma] = callback_diffeo(nu, mu{2}, Y', X{2}', rho,epsilon,options);
    
    %%% display %%%
    clf; hold on;
    % underlying grid
    opt.edge_color = [1 1 1]*.8;
    display_grid(  reshape(X0 + K(X0,Z)*a, [n0 n0 2]), opt );
    % measures
    opt.color = [1 0 0];
    plot_tensors_2d_scattered(mu{1}, X{1}, opt);
    opt.color = [0 0 1];
    plot_tensors_2d_scattered(mu{2}, X{2}, opt);  
    opt.color = [1 0 1];
    plot_tensors_2d_scattered(nu, Y, opt); 
    % add links
    E = trM(gamma,1); [~,I] = sort(E(:), 'descend'); 
    I = I(1:N(1)*2); 
    [i,j] = ind2sub(size(E), I(:)); 
    plot( [Y(i,1), X{2}(j,1)]', [Y(i,2), X{2}(j,2)]', 'color', [1 1 1]*.4 );
    drawnow;
    %%% update positions %%%
    % transfer to gradient with respect to a
    if use_a
        nabla_a = Kmat' * nabla_x';
        % additional regularization of the gradient
        if sigma_grad>0
            nabla_a = reshape(nabla_a, [p p 2]);
            nabla_a = imgaussfilt(nabla_a, sigma_grad);
            nabla_a = reshape(nabla_a, [p*p 2]);
        end
        %
        a = a - tau * nabla_a;
    else
        Y = Y - tau*nabla_x';
    end
end




