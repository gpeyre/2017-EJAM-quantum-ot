function w = sparse_marginal(i,j,gamma_ij, N, do_tranpose)

% sparse_marginal - compute row or cols sums
%
%   w = sparse_marginal(i,j,gamma_ij, N, do_tranpose)
%
%   Handles coupling computed in a sparse format as
%       gamma_ij(:,:,k) = gamma(:,:,i(k),j(k))
%   (N(1),N(2)) is the number of row/col
%
%   do_transpose=0: sums along colums
%   do_transpose=1: sums along lines
%
%   Copyright (c) 2016 Gabriel Peyre


if do_tranpose==1
    [i,j] = deal(j,i);
    N = [N(2) N(1)];
end

d = size(gamma_ij,1);
w = zeros(d,d,N(1));
for a=1:d
    for b=1:d
        H = sparse(i,j,squeeze(gamma_ij(a,b,:)), N(1), N(2) );
        w(a,b,:) = reshape(full(sum(H,2)), [1 1 N(1)]);
    end
end
    
end