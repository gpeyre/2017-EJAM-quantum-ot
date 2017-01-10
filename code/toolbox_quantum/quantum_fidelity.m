function [W,nabla_x,gamma] = quantum_fidelity(mu, nu, x, y, rho, epsilon, options)

% quantum_fidelity - compute value and gradient of Quantum OT
%
%   [W,nabla_x] = quantum_fidelity(mu,nu, x, y, options);
%
%   Compute gradient with respect to x (of size (2,N)) of W(mu,nu) where
%       mu is a tensor field at positions x
%       nu is a tensor field at positions y 
%   The cost is 1/2*|x-y|^2
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;

d = size(x,1); % embedding dimension
N = [size(x,2) size(y,2)];

% cost
c = 1/2 * tensor_id( distmat(x,y).^2, d );
% run sinkhorn
options.niter = getoptions(options, 'niter_sinkhorn', 500); % ok for .05^2
options.disp_rate = NaN;
options.tau = 1.8*epsilon/(rho+epsilon);  % prox step, use extrapolation to seed up
[gamma,u,v,err] = quantum_sinkhorn(mu,nu,c, epsilon,rho, options);
% gradient of cost
nablaC = @(x,y)repmat(x,[1 1 size(y,2)]) - ...
    repmat(reshape(y,[size(y,1) 1 size(y,2)]),[1 size(x,2)]);
% gradient with respect to positions x
nabla_x = sum( nablaC(x,y) .* ...
    repmat(reshape(trM(gamma,1), [1 N(1) N(2)]), [d 1 1]),  ...
    3 );

W = 0; % TODO.

end