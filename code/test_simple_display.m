%%
% Test for display of tensor-valued measures (does nothing but displaying ellpises).

rep = 'results/display/';
[~,~] = mkdir(rep);

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_geometry/');
addpath('toolbox_connections/');
addpath('data/meshes/');
addpath('data/images/');

%%
% Image diffusion tensors. 

n = 256;
name = 'hibiscus';
f0 = load_image(name, n); 
f0 = rescale(sum(f0,3));

op = load_helpers(n);
options.sub = 12;
clf; sigma0 = 2;
plot_tensor_field(op.rotate(op.ST(f0,sigma0)), f0, options);
saveas(gcf, [rep 'structure-tensors.png'], 'png');


%%
% Surface tensor.

name = 'moomoo';
[V,F] = read_mesh([name '.off']);
N = size(V,2);
% center/normalize the mesh
V = V - repmat(mean(V,2), [1 N]);
V = V ./ max(sqrt(sum(V.^2)));
opt.singularity = [1,2]; % (node,order), sum(order)=2
U = mesh_eigenbasis(V,F, opt);
%
clf;
opt.scale = .025;
opt.color = ones(N,1);
Id = tensor_diag(ones(N,1),.3*ones(N,1));
plot_tensor_mesh(V,F, U, Id, opt);
colormap gray(256);
saveas(gcf, [rep 'mesh.png'], 'png');

%%
% line of 2x2

name =  'aniso-iso';
N = [20 20];
mu = load_tensors_pair(name,N);
clf;
options.nb_ellipses = N(1);
plot_tensors_1d(mu{1}, options);
saveas(gcf, [rep '1d-2x2.png'], 'png');


%%
% line of 3x3

name = 'aniso-iso-3x3';
N = [20 20];
mu = load_tensors_pair(name,N);
opt.nb_ellipses = N(1);
opt.color = 1;
clf;
plot_tensors_volume(mu{1}, opt);
saveas(gcf, [rep '1d-3x3.png'], 'png');
