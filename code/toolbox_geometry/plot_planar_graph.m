function plot_planar_graph(XY,E, lw, col)

% plot_planar_graph - display planar graph with edge weight and colors.
%
% plot_planar_graph(XY,E, lw, col)
%
%   XY is a (2,N) position matrix.
%   E is (2,P) edge matrix
%   lw is (1,P) length matrix
%   col is (3,P) color matrix
%
%   Copyright (c) 2015 Gabriel Peyre

EX = reshape(XY(1,E),size(E));
EY = reshape(XY(2,E),size(E));

if size(XY,1)==2
    h = plot(EX,EY);
elseif size(XY,1)==3
    EZ = reshape(XY(3,E),size(E));
    h = plot3(EX,EY,EZ, 'k-');
else
    error('Problem');
end

for i=1:length(h)
    set(h(i), 'LineWidth', lw(i), 'color', col(:,i));
end

end