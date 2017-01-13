function [c,c0] = ground_cost(N,ground_d,tensor_d)

% ground_cost - ground cost on a uniform grid
%
%   [c,c0] = ground_cost(N,d);
%
%   Copryright (c) 2016 Gabriel Peyre

if length(N)==1
    N = [N N];
end
if nargin<3 
    tensor_d = ground_d;
end

switch ground_d
    case 1
        x = linspace(0,1,N(1));
        y = linspace(0,1,N(2));
        [Y,X] = meshgrid(y,x);
        c0 = abs(X-Y).^2;
        resh = @(x)reshape( x, [2 2 N(1) N(2)]);
        c = resh( tensor_diag(c0(:),c0(:)) );
    case {2, '2d-per'}      
        n = N(1); N = n*n;
        x = (0:n-1)'/n;
        [y,x] = meshgrid(x,x);
        [X1,X2] = meshgrid(x(:),x(:));
        [Y1,Y2] = meshgrid(y(:),y(:));
        dX = X1-X2; dY = Y1-Y2;
        if strcmp(ground_d, '2d-per')
            % periodic BC
            dX = mod(dX,1); dX(dX>1/2) = dX(dX>1/2) - 1;
            dY = mod(dY,1); dY(dY>1/2) = dY(dY>1/2) - 1;
        end
        c0 = dX.^2 + dY.^2;
        resh = @(x)reshape( x, [tensor_d tensor_d N N]);
        c = resh( tensor_id(c0(:),tensor_d) );
    otherwise 
        error('Not implemented.');
end

end