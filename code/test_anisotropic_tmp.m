% anisotropic metric
using_factorization = 1;
if not(using_factorization)
    opt.CholFactor = 0;
else
    opt.CholFactor = gamma;
end
opt.laplacian_type = 'fd'; % finite differences
if isempty(find(mask==0))
    opt.laplacian_type = 'superbases'; % J.M. Mirebeau's method
end
[blur, Delta, Grad] = blurAnisotropic(M,opt);
filtIter = 5;
K = @(x)blur(x,gamma,filtIter);

% TODO: fix naming 
% quantum_interp -> quantum_interp_single
% compute_quantum_interp -> quantum_interp