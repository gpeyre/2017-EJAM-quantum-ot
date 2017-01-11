function g = texture_synth(Ts,f0,s)

% texture_synth - gaussian texture synthesis
%  
%   g = texture_synth(Ts,f0,s)
%
%   Ts is the sqrtm of the covariance.
%
%   Copyright (c) 2017 Gabriel Peyre

if nargin<3
    s=0;
end

n = size(Ts,3);

% synthesize by convolving with a noise
if s>0
    randn('seed', s);
end
W = randn(n)/n;
W = W-mean(W(:))+1/n^2;
W = fft2(W);
%
g = tensor_mult(Ts, reshape(repmat(W, [3 1]), [3 1 n n]));
g = permute(squeeze(g), [2 3 1]);
g = real( ifft2(g) );
g = g + repmat( f0, [n n 1] );

end
