function p = periodic_comp(h)

% periodic_comp - compute the periodic component of an image
%
%   p = periodic_comp(h);
%
%   Uses the method developped by Lionel Moisan.
%
%   Copyright (c) 2016 Gabriel Peyre

[m, n, dm] = size(h);

p = zeros(size(h));
for d = 1:dm
    Hhat = intLaplace(h(:,:,d));
    FHhat = fft2(Hhat);
    
    [X Y] = meshgrid(0:n-1, 0:m-1);
    FPhat = FHhat./(4 -2.*cos(2.*X*pi/n) -2.*cos(2.*Y*pi/m));
    FPhat(1,1) = sum(sum(h(:,:,d)));
    
    p(:,:,d) = real(ifft2(FPhat));
end

end

function p = intLaplace(h)
% Compute the discrete Laplacian in the interior of the image domain
Hext = [0, h(1,:), 0; h(:,1), h, h(:,end); 0, h(end,:), 0;];

p = 4.*Hext - ( circshift(Hext,[0 1]) + circshift(Hext,[1 0]) +...
    circshift(Hext,[-1 0]) + circshift(Hext,[0 -1]) );

p = p(2:end-1, 2:end-1);

end