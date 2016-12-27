function mu = tensor_coarsening(mu,d,nsub)

% tensor_coarsening - coarsening by averaging
%
%   mu = tensor_coarsening(mu,d,nsub);
%
%   mu is reduced in size by 2^nsub along each dimension.
%   d is the number of dimensions (assumes same size along each dimension).
%
%   Copyright (c) 2016 Gabriel Peyre

if nsub<=0
    return;
end

if nsub>1
    for j=1:nsub
        mu = tensor_coarsening(mu,d,1);
    end
    return;
end

if iscell(mu)
    for k=1:length(mu)
        mu{k} = tensor_coarsening(mu{k},d,1);
    end
    return;
end

s = size(mu,1);
N = size(mu,3);
switch d
    case 2
        n = sqrt(N);
        mu = reshape( mu, [s s n n] );
        nu = zeros(s,s,n/2,n/2);
        for a=1:s
            for b=1:s
                nu(a,b,:,:) =   mu(a,b,1:2:end,1:2:end) + mu(a,b,2:2:end,2:2:end) + ...
                                mu(a,b,2:2:end,1:2:end) + mu(a,b,1:2:end,2:2:end);
            end
        end
        mu = reshape(nu, [s s n*n/4]);
    otherwise
        error('Not implemented');
end

end