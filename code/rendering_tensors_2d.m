function rendering_tensors_2d(nu,n1, naming, options)

% display of a collection of 2D tensor fields using anistotropic diffusion + ellipsoids
%
%   rendering_tensors_2d(nu,n1, naming, options);
%
%   n1 is the size of the upscaled background image.
%   Set naming=[] for no saving.
%
%   Copyright (c) Gabriel Peyre 2016

if not(iscell(nu))
    nu = {nu};
end

N = size(nu{1}, 3); 
n = floor(sqrt(N));
op = load_helpers(n);
m = length(nu);

% diffusion parameters
options.null = 0;
diffus_tau = getoptions(options, 'diffus_tau', .06); % diffusion step size
diffus_t = getoptions(options, 'diffus_t', 20); % final time

% ellipsoids parameters
opt.nb_ellipses = getoptions(options, 'nb_ellipses', 20);
opt.color = [1 0 0];
opt.color_edge = [.7 0 0];
opt.scaling = .5;

f0 = randn(n1);
% display as diffusion
for k=1:m
    % ellipsoids
    C_t = reshape(nu{k},[2 2 n n]);
    % rescale to have eigenvalues <1
    [e1,e2,l1,l2] = tensor_eigendecomp(C_t);
    C_t = C_t/max(l1(:));
    % do the diffusion
    T_t = op.C2T(C_t);
    Ti_t = image_resize(T_t,[n1 n1 3]);
    fi_t = anisotropic_diffusion(Ti_t, f0, diffus_t, diffus_tau);
    if not(isempty(naming)) % save to file
        imwrite(rescale(fi_t), [naming '-' num2str(k) '.png'], 'png');
    end
    % display with overlayed tensors
    opt.image = fi_t;
    clf;
    plot_tensors_2d(op.T2C(Ti_t), opt);
    if not(isempty(naming)) % save to file
        saveas(gcf, [naming '-ellispses-' num2str(k) '.png'], 'png');
    end
end

end