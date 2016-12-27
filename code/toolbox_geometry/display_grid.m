function display_grid(T)

% display_grid - display a warped grid
%
%   display_grid(T)
%
%   T should be of size (n,n,2)
%
%   Copyright (c) Gabriel Peyre

n = size(T,1);
surf(T(:,:,1), T(:,:,2), zeros(n));
caxis([-1 0]);
view(2); axis off; axis equal; colormap gray(256); axis ij;

end