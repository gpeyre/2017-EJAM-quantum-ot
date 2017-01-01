function f1 = anisotropic_diffusion(T, f, t, tau, options)

% anisotropic_diffusion - perform anisotropic diffusion along tensor T
%
%   f = anisotropic_diffusion(T,f,t,tau, options);
%
%   Solves : 
%       df/dt = -Delta_T(f)  s.t.  f=contrast(f)
%
%   where
%      Delta_T(f) = div(T*nabla(f))
%   is the anisotropic (Riemannian) Laplacian.
%   
%   options.laplacian='fd' select usual finite differences Laplacian, 
%       = 'superbases' selects the superbases method developped by J.M. Mirebeau
%       which gives better results for anisotropic metrics. 
%
%   T should be (n,n,3) representing entries of a 2x2 tensor field.
%
%   options.contrast is a projection function, typically histogramm
%   equalization (default: uniform histogram).
%
%   tau is step size
%   t is final time
%
%   Copyright (c) 2016 Gabriel Peyre

n = size(f,1);

options.null = 0;
verb = getoptions(options, 'verb', 1);
laplacian = getoptions(options, 'laplacian', 'fd');

% fwd differences
options.order = 1;
options.bound = 'sym';
nabla  = @(f)grad(f,options);
nablaS = @(f)-div(f,options); % adjoint

% impose contrast
x0 = reshape((1:n^2)/n^2,[n,n]);
contrast = getoptions(options, 'contrast', @(x)hist_eq(x,x0));


% Tensor x vector
Mult = @(S,u)cat(3, S(:,:,1).*u(:,:,1) + S(:,:,3).*u(:,:,2), ...
                          S(:,:,3).*u(:,:,1) + S(:,:,2).*u(:,:,2) );
 
switch laplacian
    case 'fd'
        Delta = @(f)-nablaS( Mult(T, nabla(f) ) );
    case 'superbases'
        M = zeros(n,n,2,2);
        M(:,:,1,1) = T(:,:,1);
        M(:,:,2,2) = T(:,:,2);
        M(:,:,1,2) = T(:,:,3);
        M(:,:,2,1) = T(:,:,3);
        opt.laplacian_type = 'superbases';
        [Blur, DeltaM, Grad] = blurAnisotropic(M, opt);
        Delta = @(x)reshape(-DeltaM*x(:), [n n]);
    otherwise
        error('Unknown');
end

niter=  round(t/tau);
kdisp = getoptions(options, 'disp_freq', max(2,round(niter/100)));

f1 = contrast(f);
if kdisp>0
    clf;
end
for it=1:niter
    if verb==1
        progressbar(it,niter);
    end
    f1 = f1 + tau * Delta(f1);
    f1 = contrast(f1);
    if kdisp>0 && mod(it,kdisp)==1
        imageplot(f1); drawnow;
    end
end

end