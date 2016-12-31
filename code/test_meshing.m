%%
% Test for anisotropic meshing.

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('toolbox_fast_marching/');
addpath('toolbox_quantum/tensor_logexp/');
addpath('data/images/');

name = 'circular';
name = '2d-bump-donut';
name = 'images';

rep = ['results/meshing/' name '/'];
[~,~] = mkdir(rep);

%%
% Load input measures.

n = 32; % Size for Sinkhorn computation
n1 = 256*2; % Size of the underlying grid for FM computation.

op = load_derivative_2d();
% rotate tensor by pi/2 for display
rotate = @(mu)cat(2, ...
    cat(1, mu(2,2,:,:), -mu(1,2,:,:)), ...
    cat(1, -mu(2,1,:,:), mu(1,1,:,:)) ); 


switch name
    case 'images'
        name_img = {'cartoon' 'bugsbunny'};
        Mu = {};  Mu1 = {};
        for k=1:2
            f0{k} = load_image(name_img{k}, n1);
            f0{k} = rescale(sum(f0{k},3));
            % use hessian
            sigma = 8;
            Mu0{k} = op.Hess( perform_blurring(f0{k},sigma) );
            % remapping function for the eigenvalues
            [e1,e2,l1,l2] = tensor_eigendecomp(Mu0{k});
            alpha = 1.5; beta = alpha;  % controls how much anisotropy
            vmin = 1e-4; % controls edge vs. smooth distribution
            v = max(abs(l1(:)));
            phi = @(x)vmin + abs(x/v).^alpha;
            psi = @(x)vmin + abs(x/v).^beta;
            Mu0{k} = tensor_eigenrecomp(e1,e2,phi(l1),psi(l2));
            % blur & sub-sample
            sigma = 25;
            for a=1:2
                for b=1:2
                    g = squeeze(Mu0{k}(a,b,:,:));
                    g = perform_blurring(g, sigma);
                    Mu1{k}(a,b,:,:) = g;
                    g = image_resize(g,[n n 1]);
                    Mu{k}(a,b,:,:) = reshape(g, [1 1 n n]);
                end
            end            
            mu{k} = reshape(Mu{k},[2 2 n*n]);
            % display
            clf;  
            opt.nb_ellipses = n; 
            opt.image = f0{k};
            plot_tensors_2d( rotate(Mu{k}), opt );           
            saveas(gcf, [rep 'input-images-' num2str(k) '.png'], 'png')
        end
        
    otherwise
        Mu = load_tensors_pair(name, n);
        mu = {};
        vmin = .03; % increase to avoid too rough variatons
        for k=1:2
            [e1,e2,l1,l2] = tensor_eigendecomp(Mu{k});
            Mu{k} = tensor_eigenrecomp(e1,e2,vmin+l1,vmin+l2);
            mu{k} = reshape(Mu{k},[2 2 n*n]);
        end
end

%%
% Sinkhorn.

c = ground_cost(n,2); % Ground cost
epsilon = (.08)^2;  % regularization
rho = 1;  % fidelity
options.niter = 300; % ok for .05^2
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 100; % controls number of travelling Diracs
opt.disp_tensors = 1;
nu = compute_quantum_interp(gamma, mu, m, 2, opt);

%%
% Display anisotropic mesh according to the interpolated metrics.

P = 500; % #sampling points, high density
P1 = 25; % low density
for k=1:m
    % up-sample
    op = load_helpers(n);
    T = op.C2T(reshape(nu{k}, [2 2 n n]));
    T = image_resize(T,[n1 n1 3]);
    Nu = op.T2C(T);
    % ensure positive eigenvalues    
    [e1,e2,l1,l2] = tensor_eigendecomp(Nu);
    vmin = 1e-4; 
    Nu = tensor_eigenrecomp(e1,e2,max(vmin,l1),max(vmin,l2));    
    % do the sampling
    [X,D,Vor,v] = farthesh_point(Nu,P);    
    % display    
    opt.nb_contours = 0;
    clf; disp_farthest_sampling(D,X,v, 'voronoi', opt);
    saveas(gcf, [rep 'input-voronoi-' num2str(k) '.png'], 'png');
    clf; disp_farthest_sampling(D,X,v, 'mesh');
    saveas(gcf, [rep 'input-mesh-' num2str(k) '.png'], 'png');
    % for first/last frame, overlay the background image if available
    if exist('f0') & (k==1|k==m)
        if k==1
            f = f0{1}; 
        else
            f = f0{2};
        end
        clf; 
        opt1.edge_color = [.8 .8 1];
        opt1.image = f;
        disp_farthest_sampling(D,X,v, 'mesh', opt1);
        saveas(gcf, [rep 'input-mesh-' num2str(k) '-img.png'], 'png');
    end
    % just a few samples
    if 0
    [X,D,Vor,v] = farthesh_point(Mu,P1);
    opt.nb_contours = 12;
    clf; disp_farthest_sampling(D,X,v, 'voronoi', opt);  
    saveas(gcf, [rep 'input-contours-' num2str(k) '.png'], 'png');  
    end
end

