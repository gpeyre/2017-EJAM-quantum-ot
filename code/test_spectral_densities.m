%%
% Test for matrix-valued densities.
% Need to implement a windowing estimator. 

addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('data/textures/');

name = 'grayfabric';
name = 'fabric';
name = 'stonewall';
name = 'grunge';

f = load_image(name);
f = f/255;
n0 = size(f,1);

n = 128*2;

options.estimator_type = 'window';
options.estimator_type = 'periodic';
options.samples =  200;

[synth_func,T,Ts,S,U] = texton_estimation(f, n, options);

% anisotropy ratio
A = ( S(:,:,3)-S(:,:,2) ) ./ ( S(:,:,3)+S(:,:,2) );
B = ( S(:,:,2)-S(:,:,1) ) ./ ( S(:,:,2)+S(:,:,1) );

img = @(x)imagesc(fftshift(x));
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

clf; image(synth_func(0)); 
axis image; axis off;


