%%
% Test for Sinkhorn and barycenters on 2x2 matrices on a surface.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('toolbox_geometry/');
addpath('toolbox_connections/');
addpath('data/meshes/');

name = 'moomoo';

rep = ['results/barycenters-meshes/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

%%
% Helpers.

crossp = @(a,b)[a(:,2).*b(:,3)-a(:,3).*b(:,2), -a(:,1).*b(:,3)+a(:,3).*b(:,1), a(:,1).*b(:,2)-a(:,2).*b(:,1)];
dotp3 = @(a,b)sum(a.*b,2);
normalize3 = @(x)x ./ repmat(sqrt(sum(x.^2,2)), [1 3]);
tprod1 = @(a,b)[...
    a(1,1,:).*b(1,1,:), a(1,1,:).*b(1,2,:), a(1,1,:).*b(1,3,:); ...
    a(1,2,:).*b(1,1,:), a(1,2,:).*b(1,2,:), a(1,2,:).*b(1,3,:); ...
    a(1,3,:).*b(1,1,:), a(1,3,:).*b(1,2,:), a(1,3,:).*b(1,3,:)];
tprod = @(a,b)tprod1( reshape(a', [1 3 size(a,1)]), reshape(b', [1 3 size(a,1)]) );


%%
% Load mesh.

% low res initial mesh
[V,F] = read_mesh([name '.off']);
N = size(V,2);

% center/normalize the mesh
V = V - repmat(mean(V,2), [1 N]);
V = V ./ max(sqrt(sum(V.^2)));
opt.singularity = [1,2]; % (node,order), sum(order)=2
U = mesh_eigenbasis(V,F, opt);

clf;
opt.scale = .025;
Id = tensor_diag(ones(N,1),.2*ones(N,1));
plot_tensor_mesh(V,F, U, Id, opt);
colormap gray(256)
saveas(gcf, [rep 'input-full.png'], 'png');

%%
% Geodesic distance.

if 0
    % approximate geodesic distance on the mesh using Varadhan's formula
    Mesh = load_mesh_operators(V,F);
    t_Varadhan = 1e-2;  % time step for Varadhan formula
    t_Varadhan = 1e-3;  % time step for Varadhan formula
    Z = inv( full( diag(Mesh.AreaV) + t_Varadhan * Mesh.Delta ));
    Z = -t_Varadhan*log(Z);
    D = sqrt(max(Z,0)); 
else
    % approximate using graph distance
    D = graph_distance_mesh(V,F);
end
D = D/max(D(:));

% display geodesic distance
[~,k] = min( V(2,:) ); [~,k] = max( D(k,:) ); % find cool seed point
opt.face_vertex_color = rescale( D(k,:)' );
clf;  plot_mesh(V,F, opt);
colormap parula(256);

% generate tensors, W is the normal, and V1/V2 supposed to be tangent
% C1/C2 are min/max curvature estimates
if 0
options.curvature_smoothing = 3;
[V2,V1,C2,C1,Cm,Gc,W] = compute_curvature(V,F,options);
[C1,C2] = deal(C1/max(abs(C1)), C2/max(abs(C1)));
end

%% 
% Generate two tensor fields

% find cool seed point
switch name
    case 'moomoo'
        c = [19 686];
    otherwise
        [~,c(1)] = min( V(2,:) );
        [~,c(1)] = max( D(c(1),:) );
        [~,c(2)] = max( D(c(1),:) );
end
sigma = .1;
aniso = .1;
mu = {};
for i=1:2
    a = exp( -D(:,c(i)).^2/(2*sigma^2) );
    mu{i} = tensor_diag(a(:),a(:)*aniso);
end

opt.scale = .04;
opt.nsub = 2;
for i=1:2
    clf;
    plot_tensor_mesh(V,F, U, mu{i}, opt);
    colormap parula(256);
    saveas(gcf, [rep 'input-' num2str(i) '.png'], 'png');
end

%%
% Ground cost and parallel transport.

c = reshape( tensor_diag(D(:).^2,D(:).^2), [2 2 N N]);

%%
% MaCann-like interpolation.

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

% regularization
epsilon = (.12)^2;  % large
epsilon = (.08)^2;  % medium
% marginal fidelity
rho = 1; 

options.niter = 500; % ok for .05^2
options.disp_rate = NaN;
tic;
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);
toc

% Compute interpolation using an heuristic formula.
m = 7; % number of barycenters
nu = compute_quantum_interp(gamma, mu, m, D.^2);
% rendering
opt.scale = .04;
for k=1:m
    clf;
    plot_tensor_mesh(V,F, U, nu{k}, opt);
    colormap parula(256);
    saveas(gcf, [rep 'interpol-' num2str(k) '.png'], 'png');
end

%%
% Barycenters.

% 
options.disp_rate = NaN;
options.niter = 500/2; % sinkhorn #iterates
nu = {};
for k=1:m
    t = (k-1)/(m-1);
    w = [1-t t];
    fprintf('Barycenter %d/%d:', k, m);
    [nu{k},gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options);
end
% rendering
for k=1:m
    clf;
    plot_tensor_mesh(V,F, U, nu{k}, opt);
    colormap parula(256);
    saveas(gcf, [rep 'barycenter-' num2str(k) '.png'], 'png');
end
