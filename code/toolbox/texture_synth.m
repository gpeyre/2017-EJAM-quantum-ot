function g = texture_synth(Ts,f0, U, m, s)

% texture_synth - gaussian texture synthesis
%  
%   g = texture_synth(Ts,f0, U,m, s)
%
%   Ts is the sqrtm of the covariance.
%
%   Copyright (c) 2017 Gabriel Peyre

if nargin<3
    U = eye(3);
end
if nargin<4
    m = zeros(3,1);
end
if nargin<5
    s=0;
end

d = size(Ts,1);
n = size(Ts,3);

% synthesize by convolving with a noise
if s>0
    randn('seed', s);
end
W = randn(n)/n;
W = W-mean(W(:))+1/n^2;
W = fft2(W);
%
g = tensor_mult(Ts, reshape(repmat(W, [d 1]), [d 1 n n]));
g = permute(squeeze(g), [2 3 1]);
g = real( ifft2(g) );
g = g + repmat( f0, [n n 1] );

%% add back mean/PCA %%
g = reshape( ( U(:,1:d) * reshape(g, [n*n d])' )', [n n 3] );
g = g + repmat( reshape(m,[1 1 3]), [n n 1] );


end
