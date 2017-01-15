%%
% test for DTI imaging data loading/analyzing.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('toolbox_anisotropic/');
% change this to the path containing DTI data.
addpath('/Users/gpeyre/Dropbox/work/wasserstein/wasserstein-tensor-valued/franco_brain_data');


subj = {'three' 'four'};
btype = [2000 2000];
%
subj = {'one' 'two'};
btype = [1000 2000];
%
subj = {'four' 'two'};
btype = [2000 2000];


rep = ['results/dti/' subj{1} '-' subj{2} '/'];
[~,~] = mkdir(rep);

d = 3; % run 3D code
d = 2; % run 2D code

resh = @(x)reshape(x, [size(x,1) size(x,2) size(x,3)*size(x,4)]);
iresh = @(x)reshape(x, [size(x,1) size(x,2) sqrt(size(x,3)) sqrt(size(x,3))]);
remax = @(x)x/max(x(:));

col = {};
for k=1:2
    t = k-1;
    % loading
    load(['subj_' subj{k} '_dwi_data_b' num2str(btype(k)) '_aligned_trilin_dt.mat']);
    load(['subj_' subj{k} '_anatomy']);
    N0 = [size(slice,1), size(slice,2)];
    Mu{k} = permute( reshape( slice, [N0 3 3] ), [3 4 1 2]);
    g{k} = rescale( crop(anatomySlice,N0) );

    % cropping
    n1 = 50; sub = 2;
    sx = 20 + (0:sub:n1-1);
    sy = 40 + (0:sub:n1-1);    
    n = n1/sub;
    mu{k} = resh( Mu{k}(1:d,1:d,sx,sy) );
    
    [e1,e2,l1,l2] = tensor_eigendecomp(mu{k});
    [a,e] = eigen_remaper(l1,l2,+1);
    a = hist_eq(a, linspace(0,.95,n*n));    
    e = hist_eq(e, linspace(.05,1,n*n)).^.8;
    [l1,l2] = eigen_remaper(a,e,-1);
    mu{k} = tensor_eigenrecomp(e1,e2,l1,l2);
    
    % color for display
    switch subj{k}
        case 'one'
            col{k} = [0 0 1];
            m=1:2;
        case 'two'
            col{k} = [1 0 0];
            m=2:3;
        case 'three'
            col{k} = [0 1 1];
            m = [1];
        case 'four'
            col{k} = [0 1 0];
            m=[1 3];
    end

    
    % save trace
    f{k}= g{k}; % 1-max( remax(trM(Mu{k},1)), 0);
    imwrite( remax(f{k}), [rep 'original-trace-' num2str(k) '.png'], 'png' );
    a = repmat(remax(f{k}), [1 1 3]);
    a(min(sx):max(sx),min(sy):max(sy),m) = a(min(sx):max(sx),min(sy):max(sy),m)/5;
    imwrite( remax(a), [rep 'original-trace-' num2str(k) '-roi.png'], 'png' );
    
    % display
    opt.nb_ellipses = 25;
    opt.image = trM(iresh(mu{k}),1);
    opt.color = col{k};
    clf; plot_tensors_2d(iresh(mu{k}(1:2,1:2,:,:)), opt);
    saveas(gcf, [rep 'input-' num2str(k) '.png'], 'png');
    
end


%%
% Compute the coupling using Sinkhorn. 

global logexp_fast_mode;
if d==3
    logexp_fast_mode = 1; % slow
else
    logexp_fast_mode = 4; % fast mex
end

% Ground cost
c = ground_cost(n,2,d);
% regularization
epsilon = (.08)^2;  % medium
% fidelity

rho = 1;  %medium
rho = .05; % low, usually better results

options.niter = 500; % ok for .05^2
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
fprintf('Sinkhorn: ');
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 100;
opt.disp_tensors = 1;
fprintf('Interpolating: ');
nu = quantum_interp(gamma, mu, m, 2, opt);

%%
% Display interpolation.

sx = 1:n-1; sy = 1:n-1; opt.nb_ellipses = n-1;
for k=1:m
    t = (k-1)/(m-1);
    T = iresh(nu{k}); T = T(:,:,sx,sy);
    opt.color = t*col{2} + (1-t)*col{1};
    opt.image = trM(T,1);
    clf; plot_tensors_2d(T(1:2,1:2,:,:), opt); drawnow;
    saveas(gcf, [rep 'interpol-rho' strrep(num2str(rho), '.', '') '-' num2str(k) '.png']);
end

