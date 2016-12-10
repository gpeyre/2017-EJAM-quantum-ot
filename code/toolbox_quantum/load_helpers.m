function op = load_helpers(n)

% load_helpers - load simple helper functions
%
%   op = load_helpers(n);
%
%   n is the width of a typical image.
%
%   Copyright (c) 2016 Gabriel Peyre

if nargin<1
    n = 128;
end

% centered differences
options.order = 2;
op.nabla  = @(f)grad(f,options);
op.nablaS = @(f)-div(f,options); % adjoint

op.dotp = @(x,y)sum(x(:).*y(:));
op.mynorm = @(x)norm(x(:));
% gaussian
t = [0:n/2 -n/2+1:-1];
[X2,X1] = meshgrid(t,t);
op.normalize = @(h)h/sum(h(:));
h = @(sigma)op.normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );
% convolution
op.cconv = @(f,h)real(ifft2(fft2(f).*repmat(fft2(h),[1 1 size(f,3)])));
op.blur = @(f,sigma)op.cconv(f,h(sigma));
% tensor uu'
op.tensorize = @(u)cat(3, u(:,:,1).^2, u(:,:,2).^2, u(:,:,1).*u(:,:,2));
% just for display, to allign the tensors with the edges
op.rotate = @(T)cat(3, T(:,:,2), T(:,:,1), -T(:,:,3));

% trace
op.trM = @(x)squeeze( x(1,1,:,:)+x(2,2,:,:) );


% turn into 2 x 2 x n x n format
tresh = @(x)reshape(x,[1 1 size(x,1) size(x,2)]);
op.T2C = @(T)[ tresh(T(:,:,1)), tresh(T(:,:,3)); tresh(T(:,:,3)), tresh(T(:,:,2)) ];
op.C2T = @(C)cat(3, squeeze(C(1,1,:,:)), squeeze(C(2,2,:,:)), squeeze(C(1,2,:,:)) );

% structure tensor
op.ST = @(f,sigma)op.blur( op.tensorize(op.nabla(f)), sigma);
% tensor/vector multiplication
op.tmult = @(T,u)cat(3, ...
    T(:,:,1).*u(:,:,1) + T(:,:,3).*u(:,:,2), ...
    T(:,:,3).*u(:,:,1) + T(:,:,2).*u(:,:,2) );

% structure tensor reconstruction energy
op.TensorFit.E = @(f,T0,sigma)1/2*op.mynorm(op.ST(f,sigma)-T0)^2;
op.TensorFit.nablaEi = @(H,u,sigma)2*op.nablaS( op.tmult(op.blur(H,sigma),u) );
op.TensorFit.nablaE  = @(f,t,sigma)op.TensorFit.nablaEi(op.ST(f,sigma)-t, op.nabla(f), sigma);

end