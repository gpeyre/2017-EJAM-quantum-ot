function display_grid(T, options)

% display_grid - display a warped grid
%
%   display_grid(T, opt)
%
%   T should be of size (n,n,2)
%
%   Copyright (c) Gabriel Peyre

options.null = 0;
edge_color = getoptions(options, 'edge_color', [0 0 0]);
face_alpha = getoptions(options, 'face_alpha', 0); 

n = size(T,1);
surf(T(:,:,1), T(:,:,2), zeros(n), 'EdgeColor', edge_color, 'FaceAlpha', face_alpha);
caxis([-1 0]);
view(2); axis off; axis equal; colormap gray(256); axis ij;

end