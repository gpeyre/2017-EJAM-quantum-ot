function plot_tensors_1d(Nu, options)

% plot_tensors_1d - plot tensor field along a line
%
%   plot_tensors_1d(Nu, options);
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;
q = getoptions(options, 'nb_ellipses', 10); % #ellipses
col = getoptions(options, 'color', [1 0 0]);
fill_ellipses = getoptions(options, 'fill_ellipses', 1);
y = getoptions(options, 'y', 0);

if iscell(Nu)
    clf; hold on;
    Q = length(Nu);
    for i=1:Q
        t = (i-1)/(Q-1);
        options.color = [t 0 1-t];
        options.y  = 1.8*i/q; % 2*t/q;
        plot_tensors_1d(Nu{i}, options);
    end
    return;
end

Nu = squeeze(Nu);


p = size(Nu,3);
I = round(linspace(1,p,q));
x = linspace(0,1,q);

[e1,e2,l1,l2] = tensor_eigendecomp(Nu);
theta = atan2(e2(:,:,2), e2(:,:,1));
l_max = max(l1(I));
scaling = getoptions(options, 'scaling', .8);
s = scaling*1/q;
hold on;
for k=1:q
    i = I(k);
    if strcmp(col, 'interp')
        t = (k-1)/(q-1);
        c = [t 0 1-t];
    else
        c = col;
    end
    ec = c;
    if fill_ellipses
        ellipse_fill(s*l1(i)/l_max,s*l2(i)/l_max,theta(i),x(k),y,c, ec);
    else
        ellipse(s*l1(i)/l_max,s*l2(i)/l_max,theta(i),x(k),y,c);
    end
end
axis tight; axis equal;
axis off;

end