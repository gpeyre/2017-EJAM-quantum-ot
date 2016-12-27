function z = formlet(z,c,s,a)

% formlet - parametric diffeomorphism
%
%   z = formlet(z,c,s,a)
%
%   z is an (n,n,2) matrix of (x,y) points.
%   c is the center of the warp
%   s is its width (localization radius)
%   a is its magnitude:
%       a>0: expansion, 
%       a<0: contraction
%
%   One must have a in s*[-0.1592,.1956] to obtain a valid diffeomorphism.
%
% Originally introduced in the paper:
%   On Growth and Formlets: Sparse Multi-Scale Coding of Planar Shape
%   J. H. Elder, Timothy D. Oleskiw, A. Yakubovicha, G. Peyre
%   Image and Vision Computing, 2013
%
%   Copyright (c) 2016 Gabriel Peyre

n = size(z,1);

rho = @(r,s,a)r + a*sin(2*pi*r/s).*exp(-r.^2/s^2);
norma = @(z)sqrt(sum(z.^2,3));
normalize = @(z)z ./ repmat(1e-20+norma(z), [1 1 2]);
g = @(z,s,a)normalize(z) .* repmat(rho( norma(z),s,a ), [1 1 2]);
const = @(c)cat(3,ones(n)*c(1),ones(n)*c(2));
z = const(c) +  g(z-const(c),s,a);

end