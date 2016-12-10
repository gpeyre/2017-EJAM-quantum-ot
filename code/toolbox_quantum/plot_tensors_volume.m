function plot_tensors_volume(Nu, options)

% plot_tensors_volume - plot volumetric tensor field along a line
%
%   plot_tensors_1d(Nu, options);
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;
q = getoptions(options, 'nb_ellipses', 10); % #ellipses
col = getoptions(options, 'color', 0);
Y0 = getoptions(options, 'y', 0);

if iscell(Nu)
    clf; hold on;
    Q = length(Nu);
    for i=1:Q
        t = (i-1)/(Q-1);
        options.color = t;
        options.y  = 1.8*i/q; % 2*t/q;
        plot_tensors_volume(Nu{i}, options);
    end
    return;
end

Nu = squeeze(Nu);
p = size(Nu,3);

% sub-sample
I = round(linspace(1,p,q));
Nu = Nu(:,:,I);
X0 = linspace(0,1,q);

U = []; S = [];
for i=1:q
    [U(:,:,i),s] = eig(Nu(:,:,i));
    S(:,i) = diag(s);
end
U = real(U); S = real(S);

% theta = atan2(e2(:,:,2), e2(:,:,1));
l_max = max(S(:,3));
s = .8*1/q; % scaling
r = getoptions(options, 'mesh_precision', 30); % mesh precision
hold on;
for k=1:q    
    [x,y,z] = sphere(r);    
    M = s * U(:,:,k) * [S(1,k).*x(:), S(2,k).*y(:), S(3,k).*z(:)]'; 
    x = X0(k) + reshape(M(1,:), [r+1 r+1]);
    y = Y0    + reshape(M(2,:), [r+1 r+1]);
    z = 0     + reshape(M(3,:), [r+1 r+1]);
    %
    surf(x,y,z,zeros(r+1,r+1)+col); 
end

t = linspace(0,1,256);
CM = [t;t*0;1-t]';

colormap(CM);
caxis([0 1]);
shading interp;
axis tight; axis equal;
axis off;
camlight

end