%%
% Test for deformation manipulations.
% The notation 

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_geometry/');
addpath('toolbox_quantum/tensor_logexp/');
% addpath('data/optical-flows/');

global logexp_fast_mode;
logexp_fast_mode = 1; % slow
logexp_fast_mode = 4; % fast mex

name = 'optical-flow';
name = 'localized';
rep = ['results/deformations/' name '/'];
[~,~] = mkdir(rep);

% domain width
n = 16;

%%
% Helpers.

Id = tensor_id(ones(n),2);
dotp = @(x,y)sum(x(:).*y(:));
mynorm = @(x)norm(x(:));
% div/grad, normalized to grid sampling
options.order = 2;
nabla  = @(f)n*grad(f,options);
nablaS = @(f)-n*div(f,options); % adjoint
% rotate tensor by pi/2 for display
rotate = @(mu)cat(2, ...
    cat(1, mu(2,2,:,:), -mu(1,2,:,:)), ...
    cat(1, -mu(2,1,:,:), mu(1,1,:,:)) ); 
op.rotate = @(T)cat(3, T(:,:,2), T(:,:,1), -T(:,:,3));
% Jacobian tensor of the displacement f(x)-x
%   H_{i,j} = partial_{i}(V_j)
Jac1 = @(V)cat(4, nabla(V(:,:,1)), nabla(V(:,:,2)));
Jac = @(V)permute(Jac1(V), [3 4 1 2]);
% adjoint of jacobian
Jac1S = @(W)cat(3,nablaS(W(:,:,:,1)),nablaS(W(:,:,:,2)));
% uniform grid
x = (0:n-1)'/n;
[Y,X] = meshgrid(x,x);
XY = cat(3,X,Y);


%%
% Load an input diffeomorphism pre-process optical flow.

switch name
    case 'optical-flow'
        % load/preprocess an optical flow
        filename = 'Dimetrodon';
        % load and inpaint
        V = readFlowFile(['data/optical-flows/' filename '/flow10.flo']); V(V>1e5)=Inf;
        V = diffusion_inpaint(V);
        % ROI
        w = 300; s = [50,110]; % region with discontinuity        
        w = 175; s = [2 396]; % smooth region
        % crop
        selx = s(1):s(1)+w-1; sely = s(2):s(2)+w-1;
        V = V(selx,sely,:);
        % resize
        V = image_resize(V,n,n,2);
        % load background image of first frame
        f = load_image( ['data/optical-flows/' filename '/frame10'] );
        f = rescale( mean(f,3) );
        f = f(selx,sely);
        delta = -1;
    case 'localized'
        delta = +1;     
        % first
        c = [.3 .3]; s = .5; a = s * .17;  
        T{1} = formlet(XY,c,s,a);       
        V{1} = T{1}-XY; % displacement field
        % second
        c = [.7 .7]; s = .35; a = s * -.15;  
        T{2} = formlet(XY,c,s,a);        
        V{2} = T{2}-XY; % displacement field
        % no background image
        f = [];
end


%%
% Rescale if needed to garantee diffeomorphism

for k=1:2
    H = Jac(V{k});
    if delta<0
        % Symetric Cauchy tensor and anti-symetric deviator, satisfy
        %   Id+eps*J = (Id+eps*D)(Id+eps*C) + O(eps^2)
        C = (H + permute(H, [2 1 3 4]))/2;
        [e1,e2,l1,l2] = tensor_eigendecomp(C);
        %   Increasing
        delta = -1/min(l2(:)); % limit value from linearized model
    end
    V{k} = delta*V{k};
    T{k} = XY + V{k};       
    J{k} = Id + delta*Jac(V{k});  
    J{k} = Jac(T{k});   
    % [e1,e2,l1,l2] = tensor_eigendecomp(J{k});    
    % compute polar decomposition, U is orthogonal, T is SDP. One should have
    % detU==1, otherwise non-diffeomorphic.
    [U{k},mu{k},detU{k}] = polar_svd(J{k});
    [e1,e2,l1,l2] = tensor_eigendecomp(mu{k}); 
    % for fun, remap the eigenvalues
    alpha = 1;
    vmin = .05;
    mu{k} = tensor_eigenrecomp(e1,e2,vmin + abs(l1-1).^alpha,vmin + abs(l2-1).^alpha);
end


%% 
% Displays

for k=1:2
    % display orientation of the rotation field
    a = permute( squeeze( U{k}(:,1,:,:) ), [2 3 1]);
    theta = atan2(a(:,:,2),a(:,:,1));
    clf; imagesc(theta);
    colormap jet(256);
    axis image; axis off; colorbar;
    saveas(gcf, [rep 'input-' num2str(k) '-rot.png'], 'png');
    % display grid warp
    clf; display_grid(T{k});
    saveas(gcf, [rep 'input-' num2str(k) '-grid.png'], 'png');
    % display tensors
    q = rescale(sqrt(sum(V{k}.^2,3))); % magnitude of displacement
    q = squeeze( trM(mu{k},1) );
    clf;
    opt.image = q;
    opt.nb_ellipses = min(n,32);
    plot_tensors_2d(rotate(mu{k}), opt);
    saveas(gcf, [rep 'input-' num2str(k) '-tensors.png'], 'png');
end

%%
% Transportation.

for k=1:2
    mu1{k} = reshape(mu{k},[2 2 n*n]);
    J1{k} = reshape(J{k},[2 2 n*n]);
end
% Ground cost
c = ground_cost(n,2);
% regularization
epsilon = (.06)^2;  % medium
% fidelity
rho = 1;  % medium
% prox param
options.tau = 1.8*epsilon/(epsilon+rho); 
options.niter = 1000;
options.disp_rate = NaN;
[gamma,u,v,err] = quantum_sinkhorn(mu1{1},mu1{2},c,epsilon,rho, options);
% Compute interpolation using an heuristic formula.
m = 9; % #barycenters
opt.sparse_mult = 15; % amount of diracs transfered
Ji = compute_quantum_interp(gamma, J1, m, 2, opt);

%%
% Reconstruct deformation field.

for k=1:m
    t = (k-1)/(m-1);
    Jk = permute( reshape(Ji{k}, [2 2 n n]), [3 4 1 2]);
    % solve min_h E(h)=1/2*|Jac1(h)-Jk|^2
    E = @(h)1/2*mynorm(Jac1(h)-Jk)^2;
    GradE = @(h)Jac1S(Jac1(h)-Jk);    
    % L-BFGS parameters
    flat = @(x)x(:); resh = @(x)reshape(x, [n n 2]);
    BfgsCallback = @(x)deal( E(resh(x)), flat(GradE(resh(x))) );
    opt.niter = 30;
    opt.bfgs_memory = 5;
    opt.report = @(x,v)v;
    T_init = (1-t)*T{1} + t*T{2};
    % run BFGS
    [h, R, info] = perform_bfgs(BfgsCallback, T_init(:), opt);
    Ti{k} = resh(h);
    % display as grid
    clf; display_grid(Ti{k}); drawnow;
    saveas(gcf, [rep 'interp-' num2str(k) '-grid.png'], 'png');
    % clf; display_grid(T_init); drawnow;
    % saveas(gcf, [rep 'interp-' num2str(k) '-grid.png'], 'png');
end