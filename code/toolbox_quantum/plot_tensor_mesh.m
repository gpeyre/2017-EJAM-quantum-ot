function plot_tensor_mesh(V,F, U, mu, options)

% plot_tensor_mesh - plot tensors along a mesh
%
%   plot_tensor_mesh(V,F, U, mu, options);
%
%   V(:,i) is a vertex
%   U(:,1,i) and U(:,2,i) are the tangents
%   U(:,3,i) is the normal
%   mu(:,:,i) is the 2x2 tensor expressed in the basis U(:,1:2,i)
%
%   opt.sub controls subdivision for improved display.
%
%   Copyright (c) 2016 Gabriel Peyre

N = size(V,2);

options.null = 0;
nsub = getoptions(options, 'nsub', 2);
no_mesh = getoptions(options, 'no_mesh', 0);

[E1,E2,C1,C2] = tensor_eigendecomp(mu);
[C1,C2] = deal(C1/max(C1(:)), C2/max(C1(:)));

% use as background color
f = getoptions(options, 'color', rescale(C1+C2));
color_ellipses = getoptions(options, 'color_ellipses', [1 0 0]);

if nsub==0
    Va = V; Fa = F; C1a = C1; C2a = C2; fa = f;
else
    options.sub_type = 'loop';
    options.verb = 0;
    [Va,Fa] = perform_mesh_subdivision(V, F, nsub, options);
    [fa,Fa] = perform_mesh_subdivision(f(:)', F, nsub, options);
end

% display ellpsoid
hold on;
opt.face_vertex_color = rescale(fa(:));
if no_mesh==0
    h = plot_mesh(Va,Fa, opt);
    h.FaceAlpha = 1;
end
r = getoptions(options, 'scale', .05); % scale 
offs = getoptions(options, 'offset', .005); % little offset
for i=1:N
    % principal direction of the tensor in R^3
    A = U(:,1,i) * E1(i,1) + U(:,2,i) * E1(i,2);
    B = U(:,1,i) * E2(i,1) + U(:,2,i) * E2(i,2);
    ec = 'k';
    ec = color_ellipses;
    ellipse3d_fill(r*A.*C1(i),r*B.*C2(i), V(:,i) + offs*U(:,3,i),color_ellipses,ec);
end
shading interp;
colormap jet(256);

end