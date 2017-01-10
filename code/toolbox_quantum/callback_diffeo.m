function [W,nabla_x,gamma] = callback_diffeo(nu,mu,y,x,rho,epsilon,options)

%%% displace the input measure %%%
Y = X{1} + Kmat*a;
nu = push_fwd(mu{1}, X{1}, Z, a);
%%% compute gradient with respect to positions %%%
[W,nabla_x,gamma] = quantum_fidelity(nu, mu{2}, Y', X{2}', rho, epsilon, options);

end