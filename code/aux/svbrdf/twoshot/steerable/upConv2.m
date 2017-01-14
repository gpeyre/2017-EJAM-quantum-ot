function out = upConv2(im,filt,step)

% Perform upsampling with [step] then convolution with filt
%   
% Input
%       im:     input image
%     filt:     filter
%     step:     sampling vector ([1 1] or [2 2])
%
% Output
%      out:     output image

% Upsampling
S = [size(im) 1];
im = reshape(im, [S(1) S(2) prod(S(3:end))]);

stop = [step ones(1,ndims(im)-ndims(step))] .* size(im);
tmp = zeros(stop);
tmp(1:step(1):stop(1),1:step(2):stop(2), :) = im;

% Convolution
out=imfilter(tmp,filt,'circular');% XXX must be circular?
out = reshape(out, [size(out,1) size(out,2) S(3:end)]);
