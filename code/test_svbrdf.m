%%
% test for interpolation of the gloss field for svbrdf rendering.

% put here the path where the data is
addpath('/Users/gpeyre/Dropbox/work/wasserstein/wasserstein-tensor-valued/matlab/svbrdf/data/');

addpath('svbrdf/');
addpath('svbrdf/twoshot/');
addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');

rep = 'results/svbrdf/';
[~,~] = mkdir(rep);

% low res size
n = 128; 
% high res size
n0 = 256;
% number of transported tensors
q = 2000; 

% generate viewing directions
vp = [-0.3 0.5 1]; % light vector
vp = [0.3 0.5 1]; % light vector
myview = @()twoshot_pointlight([n0 n0], vp, 1);
mycontrast = @(img)sqrt(max(0,min(1,img)));
myrender = @(mat)mycontrast( twoshot_render(mat, myview()) );

%%
% Load input fields.

names = {'leather_white' 'metal_black'};
Mat = {}; f = {}; T = {}; 
for k=1:2
    load([names{k} '.mat']);
    % crop
    Mat{k}.alb_diff = crop(mat.alb_diff, [n0 n0]);
    Mat{k}.alb_spec = crop(mat.alb_spec, [n0 n0]);
    Mat{k}.normal = crop(mat.normal, [n0 n0]);
    Mat{k}.gloss = crop(mat.gloss, [n0 n0]);
    Mat{k}.alpha = mat.alpha; 
    % down-sample
    Mat{k}.gloss1 = image_resize(Mat{k}.gloss, [n n 2 2]);
    T{k} = permute(Mat{k}.gloss1, [3 4 1 2]);
    % only retain most energetic positions
    x = (0:n-1)'/n; [Y,X] = meshgrid(x,x);
    E = trM(T{k}, 1); [~,I] = sort(E(:), 'descend'); I = I(1:q);
    xy{k} = [X(I) Y(I)];
    mu{k} = T{k}(:,:,I);
    % rendering
    f{k} = myrender(Mat{k});
    imwrite(rescale(f{k}), [rep 'input-' num2str(k), '.png'], 'png');
    % save the trace for comparison 
    a = trM(reshape(T{k}, [2 2 n n]), 1);
    imwrite(rescale(a), [rep 'input-' num2str(k), '-trace.png'], 'png');
    % display
    t = (k-1);
    opt.color = [t 0 1-t];
    opt.nb_ellipses = 16; 
    opt.image = rescale(f{k});
    clf;
    A = T{k}; 
    A = permute(Mat{k}.gloss, [3 4 1 2]);
    clf; plot_tensors_2d( A, opt);
    % re-upsample
    mat1 = Mat{k}; 
    mat1.gloss = image_resize(Mat{k}.gloss1, [n0 n0 2 2]);
    f1 = myrender(mat1);
    imwrite(rescale(f1), [rep 'input-' num2str(k), '-lr.png'], 'png');    
end

%% 
% Sinkhorn.

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

c0 = distmat(xy{1}',xy{2}');
c = reshape(tensor_diag(c0(:),c0(:)), [2 2 q q]);
epsilon = (.08)^2;  % regularization
rho = 1;  % fidelity

options.niter = 500; % ok for .05^2
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);


%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 20;
nu = quantum_interp_free(gamma, mu, xy, n, m, opt);

% compare the two
imageplot( { trM(reshape(T{1}, [2 2 n n]), 1), ...
             trM(reshape(nu{1}, [2 2 n n]), 1) });

%%
% Rendering.

% trace of input tensors
for s=1:2  
    Tr{s} = Mat{s}.gloss(:,:,1)+Mat{s}.gloss(:,:,4);
end

interpl = @(t,a,b)(1-t)*a+t*b;
normalize = @(x) x ./ repmat( sqrt(sum(x.^2,3)), [1 1 3] ); 
perm = @(T)permute(T,[3 4 1 2]);
% 
base_material = 1;
base_material = 'interp';
% gaussian smoothing to avoid artifacts
sigma = 4;
for k=1:m
    t = (k-1)/(m-1);
    switch base_material
        case 'interp'
            % use linear interpolation
            Mat1.alb_diff = (1-t)*Mat{1}.alb_diff + t*Mat{2}.alb_diff;
            Mat1.alb_spec = (1-t)*Mat{1}.alb_spec + t*Mat{2}.alb_spec;
            Mat1.normal = normalize( (1-t)*Mat{1}.normal + t*Mat{2}.normal );
        otherwise
            % use fixed interpolation
            Mat1 = Mat{base_material};
    end
    Mat1.alpha = Mat{1}.alpha;
    % save trace for comparison
    a = trM(reshape(nu{k}, [2 2 n n]), 1);
    imwrite(rescale(a), [rep 'interp-trace-' num2str(k), '.png'], 'png');
    % up-sample gloss
    G = permute( reshape( nu{k}, [2 2 n n] ), [3 4 1 2]);
    G = image_resize(G, [n0 n0 2 2]);
    % smooth
    for s=1:4
        G(:,:,s) = imgaussfilt(G(:,:,s), sigma);
    end
    % correct for possible negative tensors
    [e1,e2,l1,l2] = tensor_eigendecomp(perm(G));
    G = perm( tensor_eigenrecomp(e1,e2,max(l1,1e-2),max(l2,1e-2)) );
    % correct the trace to have a good distribution
    Tr1 = G(:,:,1)+G(:,:,4);
    a = (1-t)*sort(Tr{1}(:))+t*sort(Tr{2}(:));
    b = hist_eq(Tr1(:), a);
    for s=1:4
        G(:,:,s) = G(:,:,s).*reshape(b, [n0 n0])./Tr1;
    end
    imwrite(rescale(G(:,:,1)+G(:,:,4)), [rep 'interp-trace-' num2str(k), '.png'], 'png');
    %
    Mat1.gloss = G;
    % export
    f1 = myrender(Mat1);
    imwrite(rescale(f1), [rep 'interp-' num2str(k), '.png'], 'png');  
end