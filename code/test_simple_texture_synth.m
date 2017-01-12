%%
% Tests of gaussian texture synthesis.


addpath('toolbox/');
addpath('toolbox_quantum/');
addpath('data/textures/');

a = dir('data/textures/*.jpg');


rep = ['results/texturesynth-samples/'];
[~,~] = mkdir(rep);


perm = @(x)permute(x, [3 4 1 2]);
resh = @(x)reshape(x, [2 2 n n]);
fshift1 = @(x)[x(end/2+1:end,end/2+1:end,:,:), x(end/2+1:end,1:end/2,:,:); ...
              x(1:end/2,end/2+1:end,:,:), x(1:end/2,1:end/2,:,:)];
fshift = @(x)perm( fshift1( perm(x) ) );

remap = @(x)x.^.5;
maxdiv = @(x)x/max(x(:));

n = 256;
options.samples = 1;
for k=1:length(a)
    name = a(k).name(1:end-4);    
    f = load_image(name) / 255;   
    imwrite(f(1:n,1:n,:),[rep name '-original.png'], 'png' );
    [Ts,f0] = texton_estimation(f, n, options); 
    g = texture_synth(Ts,f0);
    imwrite(g(1:n,1:n,:),[rep name '-synthesized.png'], 'png' );
    % spectrum
    Hs = fshift( Ts );
    A = real(trM(Hs,1)); A(end/2+1,end/2+1) = Inf; A(end/2+1,end/2+1) = min(A(:));
    imwrite(maxdiv(remap(A)), [rep name '-spectrum.png'], 'png');
end