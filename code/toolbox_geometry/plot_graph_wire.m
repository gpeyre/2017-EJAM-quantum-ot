function plot_graph_wire(XY,E,xedge, xnode, options)

% plot_planar_graph - display planar graph with edge weight and colors.
%
% plot_planar_graph(XY,E, lw, col)
%
%   XY is a (2,N) position matrix.
%   E is (2,P) edge matrix
%   xedge is (1,P) signal (or []), assumed to be in [0,1]
%   xnode is (1,N) signal (or []), assumed to be in [0,1]
%
%   Copyright (c) 2015 Gabriel Peyre

options.null = 0;
edge_width = getoptions(options, 'edge_width', 5);
node_width = getoptions(options, 'node_width', 50);


if size(E,1)>size(E,2)
    E = E';
end

CM = jet(256);
CM = parula(256);

xedge = xedge(:);
xnode = xnode(:);

P = size(E,2); % #edges
N = size(XY,2); % #nodes

%% function on edges %%
if not(isempty(xedge))
    % width
    lw = .5 + edge_width*xedge;
    % colors
    c = xedge;
    Nc= size(CM,1);
    col = CM( round(c*(Nc-1))+1,:)';    
else
    lw = ones(P,1);
    col = zeros(3,P);
end
plot_planar_graph(XY,E, lw, col);

%% function on nodes %%
if not(isempty(xnode))
    % width
    ps = node_width * max(1e-10,xnode);
    % colors
    if isempty(xedge)
        c = xnode;
        Nc= size(CM,1);
        col = CM( round(c*(Nc-1))+1,:)';
    else
        col = zeros(3,N);
    end
    plot_scattered(XY, ps, col);
else
    ps = 10*ones(N,1);
	col = zeros(3,N);
end

axis tight; axis equal; axis off;
    
end

%%%%

%%%%