function plot_tensors_2d(Nu, options)

% plot_tensors_2d - plot 2D tensor field
%
%   Nu should have size (2,2,n,n)
%
%   plot_tensors_2d(Nu, options);
%
%   options.nb_ellipses controls number of ellipses along each dimension.
%   options.image displays a background image.
%
%   Copyright (c) 2016 Gabriel Peyre


options.null = 0;
q = getoptions(options, 'nb_ellipses', 10); % #ellipses
col = getoptions(options, 'color', [1 0 0]);
ecol = getoptions(options, 'color_edge', col);
fill_ellipses = getoptions(options, 'fill_ellipses', 1);
scaling = getoptions(options, 'scaling', .8);
f = getoptions(options, 'image', []);

if not(isempty(f))
   % need to exchange X/Y to correct for matlab alignement 
   Nu = permute(Nu,[1 2 4 3]);
   Nu = Nu(2:-1:1, 2:-1:1,:,:);
end

p = size(Nu,3);
I = round(linspace(1,p,q));
x = linspace(0,1,q);

[e1,e2,l1,l2] = tensor_eigendecomp(Nu);
theta = pi/2-atan2(e2(:,:,2), e2(:,:,1));
l_max = max(max(l1(I,I)));
s = scaling*1/q; % scaling
hold on;
if not(isempty(f))
    t = linspace(0,1,size(f,1));
    imagesc(t,t,f); 
    if size(f,3)==1
        colormap gray(256);
    else
        colormap jet(256);
    end
end
for k=1:q
    for l=1:q
        i = I(k); j = I(l);
        if fill_ellipses
            ellipse_fill(s*l1(i,j)/l_max,s*l2(i,j)/l_max,theta(i,j),x(k),x(l),col, ecol);
        else
            ellipse(s*l1(i,j)/l_max,s*l2(i,j)/l_max,theta(i,j),x(k),x(l),col);
        end
    end
end
axis tight; axis equal;
axis off; axis ij;

end