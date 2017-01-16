%%
% Test for Sinkhorn and barycenters on 2x2 matrices on a surface.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('toolbox_geometry/');
addpath('toolbox_connections/');
addpath('data/meshes/');

name = 'moomoo';

rep = ['results/interpolation-meshes/' name '/'];
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
T = tensor_diag(ones(N,1),.2*ones(N,1));
plot_tensor_mesh(V,F, U, T, opt);
colormap gray(256)
saveas(gcf, [rep 'input-full.png'], 'png');

%%
% Approximate geodesic distance.

D = graph_distance_mesh(V,F);
D = D/max(D(:));

% display geodesic distance
[~,k] = min( V(2,:) ); [~,k] = max( D(k,:) ); % find cool seed point
opt.face_vertex_color = rescale( D(k,:)' );
clf;  plot_mesh(V,F, opt);
colormap parula(256);

%% 
% Generate two tensor fields

% find cool seed points for the "bumps"
switch name
    case 'moomoo'
        h = {[19] [686]};
        h = {[9 374 424 524 787 804] [258]}; 
    otherwise
        [~,h(1)] = min( V(2,:) );
        [~,h(1)] = max( D(h(1),:) );
        [~,h(2)] = max( D(h(1),:) );
end

sigma = [.12 .06]; % width of each "bump"
aniso = [.03 1]; % anisotropy of each "bump"
ampl = [1 .15];
mu = {};
vmin = .02;
vmin = .06;
for k=1:2
    a = zeros(N,1);
    for i=1:length(h{k})
        a = a + exp( -D(:,h{k}(i)).^2/(2*sigma(k)^2) );
    end
    mu{k} = tensor_diag(vmin+a(:),vmin+a(:)*aniso(k));
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
% Sinkhorn.

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

% ground cost
c = reshape( tensor_diag(D(:).^2,D(:).^2), [2 2 N N]);
epsilon = (.08)^2;  % regularization
rho = 1;  % marginal fidelity
options.niter = 500; % ok for .05^2
options.disp_rate = NaN;
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

%%
%  MaCann-like interpolation.

m = 9; % number of barycenters
options.sparse_mult = 10;
nu = quantum_interp(gamma, mu, m, D.^2, options);
% correct for problems
for k=1:m
    [e1,e2,l1,l2] = tensor_eigendecomp(nu{k});
    l1 = max(l1,vmin); l2 = max(l2,vmin);
    nu{k} = tensor_eigenrecomp(e1,e2,l1,l2);    
end
    
%%
% Rendering.

opt.scale = .04;
for k=1:m
    t = (k-1)/(m-1);
    clf;
    opt.color_ellipses = [t 0 1-t];
    plot_tensor_mesh(V,F, U, nu{k}, opt);
    colormap parula(256); drawnow;
    saveas(gcf, [rep 'interpol-' num2str(k) '.png'], 'png');
end

%% 
% Anisotropic diffusion rendering. 

clear optD;
optD.niter = 150; % diffusion #iter
optD.tau = .05/1e4;
for k=1:m
    [f{k}, optD.Va,optD.Fa,T,optD.Ua] = anisotropic_diffusion_mesh(V, F,nu{k}, [], optD);
    % plot using colors
    opt.face_vertex_color = f{k};
    clf;  plot_mesh(optD.Va,optD.Fa, opt);
	colormap jet(256); drawnow;  
    saveas(gcf, [rep 'anisodiffus-' num2str(k) '.png'], 'png');
    % plot using elevations
    W = squeeze( optD.Ua(:,3,:) );
    kappa = .03; % strenght of elevation
    clf;  plot_mesh(optD.Va + kappa*repmat(f{k}'-1/2, [3 1]) .* W ,optD.Fa);
	colormap gray(256); drawnow;  
    saveas(gcf, [rep 'anisodiffus-' num2str(k) '.png'], 'png');
end
