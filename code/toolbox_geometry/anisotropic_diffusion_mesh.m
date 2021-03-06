function [f, Va,Fa,Ta,Ua] = anisotropic_diffusion_mesh(V, F, T, f, options)

% anisotropic_diffusion_mesh - perform an anisotropic tensor diffusion on a mesh
%
%   [f, Va,Fa,Ta,Ua] = anisotropic_diffusion_mesh(V, F, T, f, options);
%
%   T should be of size (2,2,N).
%
%   Copyright (c) 2017 Gabriel Peyre

options.null = 0;
Fa = getoptions(options, 'Fa', []);
Va = getoptions(options, 'Va', []);
Ta = getoptions(options, 'Ta', []);
Ua = getoptions(options, 'Ua', []);
nsub = getoptions(options, 'nsub', 2); 

% subdivide mesh 
if isempty(Va)
    options.sub_type = 'loop'; options.verb = 0;
    [Va,Fa] = perform_mesh_subdivision(V, F, nsub, options);
    opt.singularity = [1,2];
    Ua = mesh_eigenbasis(Va,Fa, opt);
end
Na = size(Va,2);
Pa = size(Fa,2);

% subdivide tensor field
if isempty(Ta)
    Ta = [];
    for i=1:2
        for j=1:2
            a = T(i,j,:);
            [a1,Fa] = perform_mesh_subdivision(a(:)', F, nsub, options);
            Ta(i,j,:) = reshape(a1, [1 1 Na]);
        end
    end
end

% ortho-bases on faces
Ua1 = ( Ua(:,:,Fa(1,:)) + Ua(:,:,Fa(2,:)) + Ua(:,:,Fa(3,:)) )/3;
Ta1 = ( Ta(:,:,Fa(1,:)) + Ta(:,:,Fa(2,:)) + Ta(:,:,Fa(3,:)) )/3;


Mesh = load_mesh_operators(Va,Fa);


contrast = getoptions(options, 'contrast', @(x)hist_eq(x,(1:Na)/Na) );
niter = getoptions(options, 'niter', 200);
tau = getoptions(options, 'tau', .2/1e4 );
kdisp = getoptions(options, 'disp_freq', max(2,round(niter/20)));

if isempty(f)
    randn('seed', 123);
    f = contrast( randn(Na,1) );
end

for it=1:niter
    g = Mesh.Grad(f);
    % coefficient in Ua basis
    Cg = tensor_mult( tensor_transp(Ua1), reshape(g, [3 1 Pa]) );
    % apply tensor
    Ch = tensor_mult( Ta1, Cg(1:2,:,:) );
    Ch(3,:,:) = 0; % 0 in normal direction
    % reconstruct in each face
    h = tensor_mult( Ua1, reshape(Ch, [3 1 Pa]) );
    u = spdiags(1./Mesh.AreaV,0,Na,Na) * Mesh.Div(h);
    f = contrast( f - tau*u );
    % display
    if kdisp>0 && mod(it,kdisp)==1
        opt.face_vertex_color = f;
        clf;  plot_mesh(Va,Fa, opt);
        colormap parula(256); drawnow;
    end
end

end
