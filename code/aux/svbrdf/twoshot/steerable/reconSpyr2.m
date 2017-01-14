function out = reconSpyr2(coeff, filter)

% Reconstruct the image from steerable pyramid
% Input:
%        coeff:     the cell structured pyramid
%       filter:     filter, typically sp3.mat, sp1.mat, sp5.mat
%
% Output:
%          out:     reconstructed image

load(filter,'lo0filt','hi0filt','lofilt','bfilts');

res = reconSpyrLevs2(coeff(2:length(coeff)),lofilt, bfilts);
temp = upConv2(res, lo0filt,[1 1]);

highpass = upConv2( coeff{1}, hi0filt,[1 1]);

out=highpass+temp;

 