%%
% Test for matrix-valued densities.
% Need to implement a windowing estimator. 

rep = 'results/texturesynthesis/';
[~,~] = mkdir(rep);

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('data/textures/');

name = 'stonewall';
name = 'grunge';

name = 'fabric';
name = 'grayfabric';
name = 'rera';
name = 'wood';
name = '5310387074_48aabd4a1b_o'; % ok
name = '8471333239_c68649960a_o'; % moyen
name = '9395669296_7d90eab03a_o'; % bof
name = 'beigefabric';

names = {'beigefabric' '5310387074_48aabd4a1b_o'};

% target size for synthesis
n = 128;
% number of transported diracs
q = 400;

img = @(x)imagesc(fftshift(x));

options.estimator_type = 'periodic';
options.samples =  200;

perm = @(x)permute(x, [3 4 1 2]);
resh = @(x)reshape(x, [3 3 n n]);

fshift1 = @(x)[x(end/2+1:end,end/2+1:end,:,:), x(end/2+1:end,1:end/2,:,:); ...
              x(1:end/2,end/2+1:end,:,:), x(1:end/2,1:end/2,:,:)];
fshift = @(x)perm( fshift1( perm(x) ) );

for k=1:2
    f{k} = load_image(names{k}) / 255;
    n0 = size(f{k},1);
    [Ts{k},f0{k},T{k},S{k},U{k}] = texton_estimation(f{k}, n, options);
    %%
    imageplot( texture_synth(Ts{k},f0{k}) );    
    clf;
    img(log(S{k}(:,:,3)+1e-10)); axis image; axis off; colorbar;    
    % only retain most energetic positions in fftshift domain.
    Hs{k} = fshift( Ts{k} );
    x = (0:n-1)'/n; [Y,X] = meshgrid(x,x);
    E = trM( Hs{k}, 1); [~,I] = sort(E(:), 'descend'); I = I(1:q);
    xy{k} = [X(I) Y(I)];
    mu{k} = Hs{k}(:,:,I);
end

%%
% Sinkhorn

% must use slow mode bc of complex matrices
global logexp_fast_mode;
logexp_fast_mode = 0;

c0 = distmat(xy{1}',xy{2}');
c = reshape(tensor_id(c0(:), 3), [3 3 q q]);
epsilon = (.08)^2;  % regularization
rho = 1;  % fidelity

options.niter = 250; 
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
[gamma,u,v,err] = quantum_sinkhorn(mu{1},mu{2},c,epsilon,rho, options);


%%
% Compute interpolation using an heuristic McCann-like formula.

m = 9;
opt.sparse_mult = 20;
nu = quantum_interp_free(gamma, mu, xy, n, m, opt);

% compare the two
imageplot( { trM(resh(T{1}), 1), ...
             trM(resh(nu{1}), 1) });
% residual
for k=1:2
    R{k} = Hs{k}-resh(nu{k});
end
% synthesize
for k=1:m
    t=(k-1)/(m-1);
    % interpolate linearly the residual
    r = (1-t)*R{1} + t*R{2};
    Hs_k = resh(nu{k}) + r;
    Rs_k = fshift(Hs_k);
    f0_k = (1-t)*f0{1} + t*f0{2};
    % display the log of energy
    A = real(trM(Hs_k,1));
    A(end/2+1,end/2+1) = Inf; A(end/2+1,end/2+1) = min(A(:));
    imwrite(rescale(log10(A+1e-5)), [rep 'spectrum-' num2str(k) '.png'], 'png');
    % synthesize the texture
    g = texture_synth(Rs_k,f0_k);
    imwrite(rescale(g), [rep 'synthesis-' num2str(k) '.png'], 'png');
end
         
return;

if 0
s = S;
% anisotropy ratio
A = ( s(:,:,3)-s(:,:,2) ) ./ ( s(:,:,3)+s(:,:,2) );
B = ( s(:,:,2)-s(:,:,1) ) ./ ( s(:,:,2)+s(:,:,1) );

clf; 
subplot(2,2,1);
img(log(S(:,:,3)+1e-10)); axis image; axis off; colorbar;
title('log(Amplitude)');
subplot(2,2,2);
img(A); axis image; axis off; colorbar;
title('Anisotropy 3/2');
subplot(2,2,3);
img(B); axis image; axis off; colorbar;
title('Anisotropy 3/2');
subplot(2,2,4);
a = U(:,1,:,:);
imageplot( squeeze(a(1,:,:,:) ) );
title('Orientation');
colormap jet(256); 

end



% clf; image(synth_func(0)); 
% axis image; axis off;


