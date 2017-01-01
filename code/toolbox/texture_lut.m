function f = texture_lut(f, name, options)

options.null = 0;

r = getoptions(options, 'cycles', 2);


saw  = @(x)mod(x*r,1);
cosa = @(x)abs(cos(pi*x*r).^2);
hat  = @(x)1-abs( 2*mod(x*r,1)-1 ); % saw-tooth



switch name
    case 'red-metal'
        remap = hat;
        c1 = [1 0 0]; c2 = [0 .5 .5];
    otherwise
        error('Unknown');
end

colorize = @(x)cat(3, ...
    (1-x)*c1(1)+x*c2(1), ...
    (1-x)*c1(2)+x*c2(2), ...
    (1-x)*c1(3)+x*c2(3) );

f = colorize( remap(f) );
        
end