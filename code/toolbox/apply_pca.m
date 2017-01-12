function [U,m,X1,c] = apply_pca(X)

% apply_pca - pca of points
%
% 	[U,m,X1,c] = apply_pca(X)
%
%   X  of size (d,n) each X(:,i) is a point.
%   U(:,k) are the principal components (leading are the first)
%   X1 = U'*X
%   c are the eigenvalues
%
%   Copyright (c) 2017 Gabriel Peyre

m = mean(X,2);
X=bsxfun(@minus,X,m);
% covariance 
C = X*X'; 
[V,c] = eig(C); % compute eignvectors
c = diag(c); % extract eigenvalues
[c,I] = sort(c, 'descend'); 
U = V(:,I); % pick leading eigenvectors
if nargout>2,
    X1=U'*X;
else
    X1=[];
end

end

