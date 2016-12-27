function W = diffusion_inpaint(V, options)

% diffusion_inpaint - fill-in missing pixels by linear diffusion
%
%    W = diffusion_inpaint(V, options);
%
%   V has size (n,p,d) and diffusion is performed along the first two
%   dimensions.
%
%   Missing values in V should be marked as Inf values.
%
%   This corresponds to solving Poisson equation inside the missing domain
%   with Dirichlet boundary conditions.
%
%   Set options.niter to be significantly larger than the width of the
%   missing region.
%
%   Copyright (c) 2016 Gabriel Peyre

average = @(x)( x([2:end end],:,:) + x([1 1:end-1],:,:) + ...
                x(:,[2:end end],:) + x(:,[1 1:end-1],:) )/4;
options.null = 0;
niter = getoptions(options, 'niter', 50);
W = V; W(isinf(W)) = 0; 
M = ~isinf(V);
for i=1:niter
    [Wold,W] = deal(W,average(W));
    W(M) = Wold(M);
end
W = W/max(W(:));

end