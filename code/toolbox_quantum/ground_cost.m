function [c,c0] = ground_cost(N,d)

% ground_cost - ground cost on a uniform grid
%
%   [c,c0] = ground_cost(N,d);
%
%   Copryright (c) 2016 Gabriel Peyre

if length(N)==1
    N = [N N];
end

switch d
    case 1
        x = linspace(0,1,N(1));
        y = linspace(0,1,N(2));
        [Y,X] = meshgrid(y,x);
        c0 = abs(X-Y).^2;
        resh = @(x)reshape( x, [2 2 N(1) N(2)]);
        c = resh( tensor_diag(c0(:),c0(:)) );
    case 2      
        n = N(1); N = n*n;
        x = (0:n-1)'/n;
        [y,x] = meshgrid(x,x);
        [X1,X2] = meshgrid(x(:),x(:));
        [Y1,Y2] = meshgrid(y(:),y(:));
        c0 = (X1-X2).^2 + (Y1-Y2).^2;
        resh = @(x)reshape( x, [2 2 N N]);
        c = resh( tensor_diag(c0(:),c0(:)) );
    otherwise 
        error('Not implemented.');
end

end