function [K,Jac,push_fwd] = load_deformation_kernel(sigma)

% load_deformation_kernel - load a RKHS and its jacobian
%
%   [K,Jac] = load_deformation_kernel(sigma)
%
%   For a set of base points Z of size (P,2), the deformation is 
%       phi_a : x -> x + K(x,Z)*a
%   where a is the coefficient vector of size P.
%
%   Jac(X,Z,a) is the jacobian of phi_a at each position defined by X of size (N,2),
%   so it is of size (2,2,N).
%
%   push_fwd(mu,X,Z,a) pushes-forward by phi_a the tensor field of mu located at 
%       X to the locations phi_a(X)
%   
%   Copyright (c) 2017 Gabriel Peyre

K = @(x,z)exp( -distmat(x',z').^2/(2*sigma^2) );
Diff = @(u,v)u*ones(1,length(v)) - ones(length(u),1) * v';
% differential with respect to x.
Kx = @(x,z)cat(3, ...
    -Diff(x(:,1),z(:,1))/(sigma^2) .* K(x,z), ...
    -Diff(x(:,2),z(:,2))/(sigma^2) .* K(x,z) );
Kz = @(x,z)-Kx(x,z);
% Jacobian of deformation field
Id = @(N)tensor_id(ones(N,1), 2);
Jac = @(X,Z,a)Id(size(X,1)) + permute( reshape( reshape( permute( Kx(X,Z), [1 3 2] ), [2*size(X,1) size(Z,1)] ) * a, [size(X,1) 2 2] ), [2 3 1]);

% Not ok
tensor_conj = @(T,D)tensor_mult( tensor_mult(D,T), tensor_transp(D) );
push_fwd1 = @(T,J)tensor_conj(T,tensor_inv(J));
push_fwd = @(T,X,Z,a)push_fwd1(T, Jac(X,Z,a) );
% Seems ok
tensor_conj = @(T,D)tensor_mult( tensor_mult(tensor_transp(D),T), (D) );
push_fwd1 = @(T,J)tensor_conj(T,(J));
push_fwd = @(T,X,Z,a)push_fwd1(T, Jac(X,Z,a) );


end