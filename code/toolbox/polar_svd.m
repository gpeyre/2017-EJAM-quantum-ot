function [U,H,detU] = polar_svd(A, isfield)
%POLAR_SVD   Canonical polar decomposition via singular value decomposition.
%   [U,H,detU] = POLAR_SVD(A) computes a matrix U of the same dimension
%   (m-by-n) as A, and a Hermitian positive semi-definite matrix H,
%   such that A = U*H.
%   U is a partial isometry with range(U^*) = range(H).
%   If A has full rank then U has orthonormal columns if m >= n
%   and orthonormal rows if m <= n.
%   U and H are computed via an SVD of A.
  
if size(A,3)>1 % modified by GP
    U = A; H = A;   
    detU = [];
    for i=1:size(A,3)
        for j=1:size(A,4)
            [U(:,:,i,j),H(:,:,i,j), detU(i,j)] = polar_svd(A(:,:,i,j));
        end
    end
    return;
end

[P,S,Q] = svd(A,'econ');
U = P*Q';
r = sum( diag(S) > norm(A,1)*eps/2 );
U = P(:,1:r)*Q(:,1:r)';
if nargout >= 2
   H = Q*S*Q';
   H = (H + H')/2;      % Force Hermitian by taking nearest Hermitian matrix.
end
if nargout >= 3
    detU = det(U);
end

end