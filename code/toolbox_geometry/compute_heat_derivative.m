function [diffDmat,D] = compute_heat_derivative(op,func,mat,x)

% compute_heat_derivative - compute adjoint of derivative for heat kernels.
%
%   diffDmat = compute_heat_derivative(op,func,mat,x);
%
%   compute 
%       D=D(x) 
%       diffDmat(u) = [partial D(x)]^*(u)
%
%   Copyright (c) 2015 Gabriel Peyre

n = size(x,1);

S = op.S(x);
W = op.W(S);
L = op.L(W);
H = op.H(L);
D = op.D(H);

diff.S = @(u)2*mat.nabla'*( (mat.nabla*x) .* repmat(u,[1 3]) );
flat = @(x)x(:);
diff.W = @(M)func.wD(S) .* flat( M( mat.E(1,:) + (mat.E(2,:)-1)*n ) );
diff.L = @(M)diag(M)*ones(1,n) - M;
diff.H = @(M)-mat.t*H*M*H;
diff.D = @(M)func.hD(H) .* M;
%
diffDmat = @(u)diff.S( diff.W( diff.L( diff.H( diff.D(u) ) ) ) );

end