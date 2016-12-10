function f1 = anisotropic_diffusion(T, f, t, tau, options)

% anisotropic_diffusion - perform anisotropic diffusion along tensor T
%
%   f = anisotropic_diffusion(T,f,t,tau, options);
%
%   Solves : 
%       df/dt = div(T*nabla(f))  s.t.  f=contrast(f)
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

   
niter=  round(t/tau);
kdisp = getoptions(options, 'disp_freq', max(2,round(niter/100)));

f1 = contrast(f);
clf;
for it=1:niter
    if verb==1
        progressbar(it,niter);
    end
    f1 = f1 - tau * nablaS( Mult(T, nabla(f1) ) );
    f1 = contrast(f1);
    if kdisp>0 && mod(it,kdisp)==1
        imageplot(f1); drawnow;
    end
end

end