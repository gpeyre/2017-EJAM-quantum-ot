function [D,Vor] = anisotropic_fm(mu,X)

% anisotropic_fm - calls the mex FM function
%
%   [D,Vor] = anisotropic_fm(mu,X);
%
%   mu should be of size (2,2,n,n)
%   X should be size (2,P)
%   D is the distance matrix
%   Vor are the vornoi cells
%
%   Copyright (c) 2016 Gabriel Peyre

n = size(mu,3);
hx = 1/n; hy = 1/n;
H = mu; 
for i=1:2
    for j=1:2
        H(:,:,i,j) = [1,0;0,1] * 10000000;
    end
end
H = permute(mu, [3 4 1 2]);
[D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, X);
Vor = Vor+1;
end