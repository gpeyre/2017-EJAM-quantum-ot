%%
% Test for matrix-valued densities.
% Need to implement a windowing estimator. 


addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('data/textures/');


names = {'beigefabric' 'fabric1'};
names = {'fabric-mix' 'stonewall'};
names = {'wood' 'fabric'};
names = {'rera' 'grunge'};
names = {'rera' 'fabric'};

rep = ['results/texturesynth/' names{1} '-' names{2} '/'];
[~,~] = mkdir(rep);

% target size for synthesis
n = 128*2;
% number of transported diracs
q = 400*3;

img = @(x)imagesc(fftshift(x));

options.estimator_type = 'periodic';
% number of samples to populate the covariance, increase to avoid rank
% defficiency problems.
options.samples =  1;

perm = @(x)permute(x, [3 4 1 2]);
resh = @(x)reshape(x, [2 2 n n]);

fshift1 = @(x)[x(end/2+1:end,end/2+1:end,:,:), x(end/2+1:end,1:end/2,:,:); ...
              x(1:end/2,end/2+1:end,:,:), x(1:end/2,1:end/2,:,:)];
fshift = @(x)perm( fshift1( perm(x) ) );

for k=1:2
    f{k} = load_image(names{k}) / 255;
    % full 3D texton
    [TsF{k},f0F{k}] = texton_estimation(f{k}, n, options);
    % dimension reduction to 2 colors    
    s = size(f{k});
    [U{k},m0{k},X,c] = apply_pca( reshape(f{k}, [prod(s(1:2)) 3])' );
    f1{k} = reshape(X(1:2,:)', [s(1:2) 2]); 
    % save input and PCA reduced input
    imwrite(clamp(f{k}(1:n,1:n,:)),[rep 'original-' num2str(k) '.png'], 'png' );
    g{k} = reshape( ( U{k}(:,1:2) * reshape(f1{k}(1:n,1:n,:), [n*n 2])' )', [n n 3] );
    g{k} = g{k} + repmat( reshape(m0{k},[1 1 3]), [n n 1] );
    imwrite(clamp(g{k}),[rep 'original-' num2str(k) '-pca.png'], 'png' );
    % estimate texton
    [Ts{k},f0{k},T{k},S{k}] = texton_estimation(f1{k}, n, options);
    Hs{k} = fshift( Ts{k} );
    % display synthesis
    clf; imageplot( texture_synth(Ts{k},f0{k}, U{k}, m0{k}) );  
    % display spectrum
    A = real(trM(Hs{k},1)); A(end/2+1,end/2+1) = Inf; A(end/2+1,end/2+1) = min(A(:));
    clf; imagesc(log(A+1e-10)); axis image; axis off; colorbar;    
    imwrite(rescale(log10(A+1e-5)), [rep 'spectrum-input-' num2str(k) '.png'], 'png');
    % only retain most energetic positions in fftshift domain.    
    x = (0:n-1)'/n; [Y,X] = meshgrid(x,x);
    E = trM( Hs{k}, 1); [~,I] = sort(E(:), 'descend'); I = I(1:q);
    xy{k} = [X(I) Y(I)];
    mu{k} = Hs{k}(:,:,I);
    % avoid singularity
    mu{k} = mu{k} + 1e-3*max(A(:)) * tensor_id(ones(q,1), 2);
end


%%
% Sinkhorn

global logexp_fast_mode;
logexp_fast_mode = 1;

c0 = distmat(xy{1}',xy{2}');
c = reshape(tensor_id(c0(:), 2), [2 2 q q]);
epsilon = (.08)^2;  % regularization
rho = 1;  % fidelity

options.niter = 250; 
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
fprintf('Sinkhorn: ');
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);

%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 20;
fprintf('Interpolation: ');
nu = quantum_interp_free(gamma, mu, xy, n, m, opt);

% residual
R = { Hs{1}-resh(nu{1}), Hs{end}-resh(nu{end}) };
% synthesize
lint = @(t,a,b)(1-t)*a + t*b;
remap = @(x)log10(A+1e-5);
remap = @(x)x.^.5;
maxdiv = @(x)x/max(x(:));
for k=1:m
    t=(k-1)/(m-1);
    % interpolate linearly the residual
    r = (1-t)*R{1} + t*R{2};
    Hs_k = resh(nu{k}) + r;
    Ts_k = fshift(Hs_k);
    f0_k = lint(t, f0{1}, f0{2} );
    m_k = lint(t, m0{1}, m0{2} );
    U_k = orthogonalize_mat( lint(t, U{1}, U{2} ) );
    % display the log of energy
    A = real(trM(Hs_k,1)); A(end/2+1,end/2+1) = Inf; A(end/2+1,end/2+1) = min(A(:));
    imwrite(maxdiv(remap(A)), [rep 'spectrum-' num2str(k) '.png'], 'png');
    % synthesize the texture
    f_k = texture_synth(Ts_k,f0_k,U_k,m_k, 123);
    imwrite(clamp(f_k), [rep 'synthesis-' num2str(k) '.png'], 'png');
end
 

imageplot( {f_k texture_synth(Ts{2},f0{2},U{2},m0{2}, 123), ...
        texture_synth(TsF{2},f0F{2},eye(3),zeros(3,1), 123) } );