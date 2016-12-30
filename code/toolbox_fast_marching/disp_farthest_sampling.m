function disp_farthest_sampling(D,X,v, disp_mode, options)

% disp_farthest_sampling - display sampling
%
%   disp_farthest_sampling(D,X,v, disp_mode);
%
%   Set disp_mode to either 'voronoi' or 'mesh'
%
%   Copyright (c) 2016 Gabriel Peyre

if nargin<4
    disp_mode = 'mesh';
end
n = size(D,1);

options.null = 0;
nb_contours = getoptions(options, 'nb_contours', 14);
f = getoptions(options, 'image', []);

hold on;
switch disp_mode
    case 'voronoi'
        % display Vornoi cells
        imagesc(1:n,1:n,D); colormap jet(256);
        if nb_contours>0
            contour(1:n,1:n,D, nb_contours, 'k');
        end
        plot(X(2,:), X(1,:), 'r.', 'MarkerSize', 15);
        axis ij; axis image; axis off; drawnow;
    case 'mesh'
        edge_color = getoptions(options, 'edge_color', 'k');
        lw = getoptions(options, 'lw', 1);
        % display mesh
        if not(isempty(f))
            imagesc(1:n,1:n,f); colormap gray(256);
        end        
        h = plot(   ...
            [X(2,v(:,1));X(2,v(:,2))], ...
            [X(1,v(:,1));X(1,v(:,2))], ...            
            'color',  edge_color, 'LineWidth', lw);
        plot([1 n n 1 1], [1 1 n n 1], 'color',  edge_color, 'LineWidth', lw);
        plot(X(2,:), X(1,:),  'r.', 'MarkerSize', 15);
        axis ij; axis image; axis off; drawnow;
end

end