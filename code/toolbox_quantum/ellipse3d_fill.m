function ellipse3d_fill(a,b,m,C,EC)

% ellipse3d_fill - draw filled 2D ellipse in 2D space
%
%   ellipse3d_fill(a,b,m,C,EC)
%
%   m in R^2 is the center
%   a,b in R^2 are the axis.
%
%   Copyright (c) 2016 Gabriel Peyre

if nargin<4
    C = 'k';
end
if nargin<5
    EC = 'k';
end

q = 32;
t = linspace(0,2*pi,q);

x = m(:)*ones(1,q) + a(:)*cos(t)+b(:)*sin(t);

if not(isempty(EC))
%    patch(x(1,:),x(2,:),x(3,:), C, 'EdgeColor', EC);
    patch('XData', x(1,:),'YData', x(2,:),'ZData', x(3,:), 'FaceVertexCData', (C(:)*ones(1,q))', 'EdgeColor', EC);    
else
    patch(x(1,:),x(2,:),x(3,:),C);
end    
end