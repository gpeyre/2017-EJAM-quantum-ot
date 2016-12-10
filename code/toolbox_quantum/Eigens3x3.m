function [Eigval,Eigvect] = Eigens3x3(A)

% Eigens3x3 - computes eigenvalues and eigenvectors of a collection of 3x3 matrices.
%
%   [Va,Ve] = Eigens3x3(A)
%
% A should be of size (3,3,N) and A(:,:,i) must be symmetric.
%
%   Copyright (c) 2016 Gabriel Peyre

crossp = @(a,b)[a(2,1,:).*b(3,1,:)-a(3,1,:).*b(2,1,:); ...
    -a(1,1,:).*b(3,1,:)+a(3,1,:).*b(1,1,:); ...
    a(1,1,:).*b(2,1,:)-a(2,1,:).*b(1,1,:)];

N = size(A,3);

% trace
trA = A(1,1,:)+A(2,2,:)+A(3,3,:);
% det
detA = A(1,1,:).*A(2,2,:).*A(3,3,:)+...
       A(1,3,:).*A(2,1,:).*A(3,2,:)+...
       A(3,1,:).*A(1,2,:).*A(2,3,:)-...
       A(1,3,:).*A(2,2,:).*A(3,1,:)-...
       A(1,1,:).*A(3,2,:).*A(2,3,:)-...
       A(3,3,:).*A(1,2,:).*A(2,1,:);
tres = -A(2,2,:).*A(3,3,:)-A(2,2,:).*A(1,1,:)-A(1,1,:).*A(3,3,:)+A(1,3,:).*A(3,1,:)+A(2,3,:).*A(3,2,:)+A(1,2,:).*A(2,1,:);

% characteristic polynomial
P = reshape([ -ones(1,1,N) trA tres detA ],4,N);
P(abs(double(P))<=eps) = 0;
% roots
Eigval = reshape(poly_root(P),3,1,N);
Eigval = real(Eigval); %eig should be real

% first eigenvector : [(A-v1*Id)*e1] ^ [(A-v1*Id)*e2]
a = A(:,1,:); a(1,1,:) = a(1,1,:) - Eigval(1,1,:);
b = A(:,2,:); b(2,1,:) = b(2,1,:) - Eigval(1,1,:);
u = crossp(a,b);
% second eigenvector : [(A-v2*Id)*e1] ^ [(A-v2*Id)*e2]
a = A(:,1,:); a(1,1,:) = a(1,1,:) - Eigval(2,1,:);
b = A(:,2,:); b(2,1,:) = b(2,1,:) - Eigval(2,1,:);
v = crossp(a,b);
% third eigenvector : u ^ v
w = crossp(u,v);

% normalize
Nu = sqrt(sum(u.^2));
Nv = sqrt(sum(v.^2));
Nw = sqrt(sum(w.^2));
Eigvect = zeros(3,3,N);
Eigvect(:,1,:) = u ./ repmat( Nu, [3 1 1] );
Eigvect(:,2,:) = v ./ repmat( Nv, [3 1 1] );
Eigvect(:,3,:) = w ./ repmat( Nw, [3 1 1] );
% detect problematic cases, hand fixing
tol = 1e-4;
I = find( (Nu<tol) | (Nv<tol) );
for i=I(:)'
    [U,D] = eig(A(:,:,i));
    Eigvect(:,:,i) = U;
    Eigval(:,:,i) = diag(D);
end

end
