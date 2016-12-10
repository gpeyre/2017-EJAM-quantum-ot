function [nu,gamma,err] = quantum_barycenters_hard(mu,c,rho,epsilon,w,options)

% quantum_barycenters_hard - Quantum barycenter computations
%
%   [nu,gamma,err] = quantum_barycenters(mu,c,rho,epsilon,w,options)
%
%   Solve for 
%       min_nu sum_i w_i W_{epsilon,c{i}}(mu{i},nu)
%   where c{i} is a cost matrix between mu{i} (rows) and nu (cols). 
%
%   options.niter is the number of Sinkhorn iterates.
%   options.niter_fixed is the number of inner loop fixed point iterates. 
%   options.display_func(gamma,mu,nu) is a callback for display.
%
%   Copyright (c) 2016 Gabriel Peyre


options.null = 0;
niter = getoptions(options, 'niter', 100);

% display marginals
disp_func = getoptions(options, 'disp_func', @mydisplay);
disp_rate = getoptions(options, 'disp_rate', max(2,ceil(niter/100)) );
expMfunc = getoptions(options, 'expMfunc', @expM); options.expMfunc = expMfunc;
logMfunc = getoptions(options, 'logMfunc', @logM); options.logMfunc = logMfunc;

% number of inputs
L = length(w);

% dimension
d = size(mu{1},1);

% size of barycenter
N0 = size(c,4);
N = [N0;N0];

% prox param
lambda = rho/epsilon;
nu = mu{1}; % ones(2,2,N0);

% relaxation parameters
tau = getoptions(options, 'tau', 1/(lambda+1));
over_iterations = getoptions(options, 'over_iterations', 1);

% init
for i=1:L
    v{i} = nu*0;
    u{i} = mu{i}*0;
    Lmu{i} = logMfunc(mu{i}); % NOT NEEDED
end
Lnu = logMfunc(nu); % NOT NEEDED


K = @(u,v) -c/epsilon ...
    -    repmat( reshape(v,[d d 1 N(2)]),[1 1 N(1) 1])/epsilon ...
    -rho*repmat( u,[1 1 1 N(2)])/epsilon;
perm = @(x)permute(x, [1 2 4 3]); 
relax = @(x,y,tau)x*(1-tau)+y*tau;
mynorm = @(x)norm(x(:));

err = [];
for it=1:niter    
    progressbar(it,niter);
    %%% update u %%%
    err(it,1) = 0;
    for i=1:L
        u1 = relax(u{i}, lse(K(u{i},v{i})) - Lmu{i}, tau);
        err(it,1) = err(it,1) + mynorm(u{i}-u1);
        u{i} = u1;
    end
    % display   
    if not(isempty(disp_func)) && mod(it,disp_rate)==1
        gamma = recomp_coupling(c,u,v,epsilon,lambda, options);
        nu = expMfunc(Lnu);
        disp_func(gamma,mu,nu);
    end

    Lnu1 = Lnu; v1 = v;    
    for it1=1:over_iterations
        %%% Stores the LSE %%%
        for i=1:L
            LSEv{i} = lse(perm(K(u{i},v1{i})));
        end        
        %%% update Lnu %%%
        Lnu1 = Lnu1*0;
        for i=1:L
            Lnu1 = Lnu1 + w(i)*( LSEv{i} + v1{i}/epsilon );
        end      
        %%% update v %%%
        for i=1:L
            v1{i} = epsilon * LSEv{i} + v1{i} - epsilon * Lnu1;
        end
    end   
    
    err(it,2) = 0;
    for i=1:L
        err(it,2) = err(it,2) + mynorm(v{i}-v1{i});
    end
    err(it,3) = mynorm(Lnu-Lnu1);
    Lnu = Lnu1;
    v = v1;  
    
    % display    
    if not(isempty(disp_func)) && mod(it,disp_rate)==ceil(disp_rate/2)
        gamma = recomp_coupling(c,u,v,epsilon,lambda, options);
        nu = expMfunc(Lnu);
        disp_func(gamma,mu,nu);
    end
end

nu = expMfunc(Lnu);
gamma = recomp_coupling(c,u,v,epsilon,lambda, options);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mydisplay(gamma,mu,nu)

L = length(gamma);
trM1 = @(x)squeeze(trM(x));

% compares trace
clf;
for i=1:L
    A0{i} = trM1(mu{i});
    B0 = trM1(nu);
    A{i} = trM1(sum(gamma{i},4));
    B{i} =  trM1(sum(gamma{i},3));
    subplot(L,2,(i-1)*L + 1);
    plot([A{i} A0{i}]); legend(['gamma\{' num2str(i) '\}_1'], ['mu\{' num2str(i) '\}']);
    axis tight;
    subplot(L,2,(i-1)*L + 2);
    plot([B{i} B0]); legend(['gamma\{' num2str(i) '\}_2'], 'nu');
    axis tight;
end
drawnow;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma = recomp_coupling(c,u,v,epsilon,lambda, options)

d = size(c,1); % dimension

if iscell(u)
    for i=1:length(u)
        gamma{i} = recomp_coupling(c,u{i},v{i},epsilon,lambda, options);
    end
    return;
end

N(1) = size(c,3); N(2) = size(c,4);
resh = @(x)reshape(x, [d d N(1) N(2)]);
flat = @(x)reshape(x, [d d N(1)*N(2)]);
a = -c/epsilon -lambda*repmat( reshape(v,[d d 1 N(2)]),[1 1 N(1) 1]) ...
    - lambda*repmat( u,[1 1 1 N(2)]);
gamma =  resh( options.expMfunc(flat(a))  );

end