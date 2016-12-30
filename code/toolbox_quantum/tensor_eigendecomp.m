function [e1,e2,l1,l2] = tensor_eigendecomp(T)

% tensor_eigendecomp - compute eigenvectors and eigenvalues of 2x2 matrices.
%
%   [e1,e2,l1,l2] = tensor_eigendecomp(T);
%
%   T is of size 2x2xNxP
%   e1 and e2 are of size NxPx2
%   l1 and l2 are of size NxP
%   eigenvalues are ordered so that l1>l2
%
%   Copyright (c) 2016 Gabriel Peyre

K11 = squeeze(T(1,1,:,:));
K12 = squeeze(T(1,2,:,:));
K21 = squeeze(T(2,1,:,:));
K22 = squeeze(T(2,2,:,:));

[n,p] = size(K11);

e1 = zeros(n,p,2);
e2 = zeros(n,p,2);
l1 = zeros(n,p);
l2 = zeros(n,p);

% trace/2
t = (K11+K22)/2;

a = K11 - t;
b = K12;

ab2 = sqrt(a.^2+b.^2);
l1 = ab2  + t;
l2 = -ab2 + t;

theta = atan2( ab2-a, b );

e1(:,:,1) = cos(theta);
e1(:,:,2) = sin(theta);
e2(:,:,1) = -sin(theta);
e2(:,:,2) = cos(theta);

end