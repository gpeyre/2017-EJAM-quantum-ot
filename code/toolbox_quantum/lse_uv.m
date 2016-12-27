function w = lse_uv(c,u,v,epsilon,rho, do_transpose)

% lse_uv - implements the LSE operator from (u,v)
%
%   w = lse_uv(c,u,v,epsilon,rho, do_transpose);
%
%   Can handle a sparse cost c.
%   do_transpose=0: sums along colums
%   do_transpose=1: sums along lines
%
%   Copyright (c) 2016 Gabriel Peyre

d = size(u,1);
N = [size(u,3), size(v,3)];
lambda = rho/epsilon;

if issparse(c)
    % -c_ij/epsilon-lambda*u_i-lambda*v_i
    [i,j,cij] = find(c);
    M = expM( -tensor_id(cij,d)/epsilon-lambda*u(:,:,i)-lambda*v(:,:,j) );
    w = sparse_marginal(i,j, M, N, do_transpose);   
    w = logM(w);
else
    perm = @(x)permute(x, [1 2 4 3]); 
    K = @(u,v) -c/epsilon ...
        -lambda*repmat( reshape(v,[d d 1 N(2)]),[1 1 N(1) 1]) ...
        -lambda*repmat( u,[1 1 1 N(2)]);
    if do_transpose==0
        w = lse(K(u,v));
    else
        w = lse(perm(K(u,v)));
    end
        
end
