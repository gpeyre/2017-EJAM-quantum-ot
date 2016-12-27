%%
% Test for Sinkhorn and barycenters on 2x2 matrices on a surface.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('toolbox_geometry/');
addpath('toolbox_connections/');
addpath('data/meshes/');

name = 'moomoo'; % 1045
name = 'lioness'; % 3400 // not spherical
name = 'gorilla'; % 2044
name = 'tre_twist'; % 800 % 1D
name = 'bull'; % 502 view(5,10)
name = 'horse'; % 417


rep = ['results/barycenter-meshes/' name '/'];
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
if size(V,2)<500
    options.sub_type = 'loop';
    options.verb = 0;
    [V,F] = perform_mesh_subdivision(V, F, 1, options);
end
N = size(V,2);

clf;  plot_mesh(V,F);

% center/normalize the mesh
V = V - repmat(mean(V,2), [1 N]);
V = V ./ max(sqrt(sum(V.^2)));
opt.singularity = [1,2]; % (node,order), sum(order)=2
U = mesh_eigenbasis(V,F, opt);

clf;
opt.scale = .025;
opt.nsub = 1; % 2 for high res
opt.offset = .002;
Id = tensor_diag(ones(N,1),.2*ones(N,1));
opt.color = ones(N,1);
opt.color_ellipses = [1 0 0];
plot_tensor_mesh(V,F, U, Id, opt);
colormap gray(256);
saveas(gcf, [rep 'input-full.png'], 'png');

%%
% Geodesic distance.

% approximate using graph distance
D = graph_distance_mesh(V,F);
D = D/max(D(:));

% display geodesic distance
[~,k] = min( V(2,:) ); [~,k] = max( D(k,:) ); % find cool seed point
opt.face_vertex_color = rescale( D(k,:)' );
clf;  plot_mesh(V,F, opt);
colormap parula(256);

if 0
clf;  plot_mesh(V,F);
idx = selectPoint(V)
end

%% 
% Generate two tensor fields

% find cool seed point
switch name
    case 'moomoo'
        c = [21 281 757 85];
    case 'horse'
        c = [207 1147 671 1217];
    otherwise
        error('Unknown');
end
sigma = .07;
aniso = .2;
mu = {};
for i=1:length(c)
    a = exp( -D(:,c(i)).^2/(2*sigma^2) );
    mu{i} = tensor_diag(a(:),a(:)*aniso);
end

opt.scale = .04;
col = [1 0 0; 0 1 0; 0 0 1; .8 .8 0];
for i=1:length(mu)
    opt.color_ellipses = col(i,:);
    opt.color = ones(N,1);
    clf;
    plot_tensor_mesh(V,F, U, mu{i}, opt);
    colormap gray(256); drawnow;
    saveas(gcf, [rep 'input-' num2str(i) '.png'], 'png');
end

opt.scale = .04;
col = [1 0 0; 0 1 0; 0 0 1; .8 .8 0];
clf; hold on;
for i=1:length(mu)
    opt.color_ellipses = col(i,:);
    opt.color = ones(N,1);
    opt.no_mesh = 0;
    if i>1
        opt.no_mesh = 1;
    end
    plot_tensor_mesh(V,F, U, mu{i}, opt);
    colormap gray(256); drawnow;
end
opt.no_mesh = 0;
saveas(gcf, [rep 'input-all.png'], 'png');



%%
% Ground cost and parallel transport.

c = reshape( tensor_diag(D(:).^2,D(:).^2), [2 2 N N]);

%%
% MaCann-like interpolation.


% regularization
epsilon = (.12)^2;  % large
epsilon = (.08)^2;  % medium
% marginal fidelity
rho = 1; 


%%
% Barycenters.

options.disp_rate = NaN;
options.niter = 80; % sinkhorn #iterates
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
        opt.color_ellipses = sum( col .* repmat(w', [1 3]) );
        clf;
        plot_tensor_mesh(V,F, U, nu{k1,k2}, opt);
        colormap gray(256); drawnow;
        saveas(gcf, [rep 'barycenter-' num2str(k1) '-' num2str(k2) '.png'], 'png');
    end
end

