function f1 = struct_tensor_reconstruc(f_init, T, sigma, options)

% struct_tensor_reconstruc - peform reconstruction from tensor field
%
%   f1 = struct_tensor_reconstruc(f_init, T, sigma, options);
%
%   Run BFGS on the energy 
%       min_f 1/2*|ST(f)-T|^2 + lambda*TV(f)
%   where lambda (the TV penalty weight) is specified in options.lambda. 
%
%   Copyright (c) 2016 Gabriel Peyre

n0 = size(f_init,1);
op = load_helpers(n0);


% TV regularizer
kappa = 0.4*1e-2;
lambda = getoptions(options, 'lambda', 0.06);
% div/grad
options.order = 1;
nabla  = @(f)grad(f,options);
nablaS = @(f)-div(f,options); % adjoint
% smoothed TV norm and its gradient
TVn  = @(f)sqrt( kappa^2 + sum(nabla(f).^2,3) );
TV = @(f)sum(sum(TVn(f)));
TVg = @(f)nablaS( nabla(f)./repmat(TVn(f), [1 1 2])  );


% L-BFGS parameters
flat = @(x)x(:); resh = @(x)reshape(x, [n0 n0]);
BfgsCallback = @(x,t,sigma)deal( op.TensorFit.E(resh(x),t,sigma) + lambda*TV(resh(x)), ...
    flat(op.TensorFit.nablaE(resh(x),t,sigma) + lambda*TVg(resh(x)) ) );
opt.niter = getoptions(options, 'niter_bfgs', 30);
opt.bfgs_memory = 5;
opt.report = @(x,v)v;


if 0
    % test on denoising
    f_init = f0{1};
    TVdenoise = @(f)deal(1/2*norm(f-f_init, 'fro')^2 + lambda * TV(f), ...
        flat( f-f_init + lambda*TVg(f) ) );
    [f, R, info] = perform_bfgs(@(x)TVdenoise(resh(x)), f_init(:), opt);
    f = resh(f);
end

% run the BFGS
[f, R, info] = perform_bfgs(@(x)BfgsCallback(x,T,sigma), f_init(:), opt);
f1 = resh(f);

end