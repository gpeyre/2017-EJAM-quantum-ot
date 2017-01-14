function coeff = buildSCFpyr2(im, ht)
% This is a modified version of buildSFpyr, that constructs a
% complex-valued steerable pyramid  using Hilbert-transform pairs
% of filters.  Note that the imaginary parts will *not* be steerable.
%
% To reconstruct from this representation, either call reconSFpyr
% on the real part of the pyramid, *or* call reconSCFpyr which will
% use both real and imaginary parts (forcing analyticity).
% Original code: Eero Simoncelli, 5/97.
% Modified by Javier Portilla to return complex (quadrature pair) channels,
% 9/97.


nbands = 4;

%============================Create lo0 and hi0 mask =================================
S = [size(im) 1];
m = S(1);
n = S(2);
%[m,n]=size(im);
[log_rad,angle]=base(m,n);

% Radial transition function (a raised cosine in log-frequency):
[Xrcos,Yrcos] = rcosFn(1,-1/2);
Yrcos = sqrt(Yrcos);
YIrcos = sqrt(1.0 - Yrcos.^2);

lo0mask = pointOp(log_rad, YIrcos, Xrcos);
hi0mask = pointOp(log_rad, Yrcos, Xrcos);
lo0mask = repmat(lo0mask, [1 1 S(3:end)]);
hi0mask = repmat(hi0mask, [1 1 S(3:end)]);
% ========================= Recursive functions ================================
imdft = fftshift(fft2(im));
lo0dft =  imdft .* lo0mask;
coeff = buildSCFpyrLevs2(lo0dft, log_rad, Xrcos, Yrcos, angle, ht-1, nbands);

% ========================= High pass residue ==============================
hi0dft =  imdft .* hi0mask;
hi0 = ifft2(ifftshift(hi0dft));

% ========================== Put together in the cell output =============
coeff = [real(hi0) coeff];