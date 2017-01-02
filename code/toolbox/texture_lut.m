function g = texture_lut(f, name, options)

% texture_lut - apply a look up table to a grayscale image to render it in color
%
%   g = texture_lut(f, name, options);
%
%   Supported name are 
%       'none' (grayscale)
%       'red-metal'
%       'wood'
%       'rust'
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;

r = getoptions(options, 'cycles', 2);

saw  = @(x)mod(x,1);
cosa = @(x)1-abs(cos(pi*x).^2);
hat  = @(x)1-abs( 2*mod(x,1)-1 ); % saw-tooth


colorize = @(x,c1,c2)cat(3, ...
    (1-x)*c1(1)+x*c2(1), ...
    (1-x)*c1(2)+x*c2(2), ...
    (1-x)*c1(3)+x*c2(3) );

switch name
    case 'none'
        g = repmat(f, [1 1 3]);
    case 'red-metal'
        c1 = [1 0 0]; c2 = [0 .3 .3];
        g = colorize( hat(f*1), c1,c2 );
    case 'wood'
        c1 = [169 100 76]/255; c2 = [45 33 30]/255;
        g = colorize( saw(f*2), c1,c2 );
    case 'rust'
        q = .7; % cutoff
        c1 = [.5 .5 .5]; c2 = [.9 .9 .9];
        g1 = colorize( f/q, c1,c2 );
        c1 = [113 64 45]/255; c2 = [83 54 19]/255;
        g2 = colorize( (f-q)/(1-q), c1,c2 );
        I = find(f<q); J = find(f>q); 
        n = size(f,1);
        g = zeros([n n 3]);
        g([I,I+n*n,I+2*n*n]) = g1([I,I+n*n,I+2*n*n]);
        g([J,J+n*n,J+2*n*n]) = g2([J,J+n*n,J+2*n*n]);
    otherwise
        error('Unknown');
end
        
end