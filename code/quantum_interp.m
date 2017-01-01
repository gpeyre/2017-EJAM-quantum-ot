function nu = quantum_interp(gamma, mu, m, d, options)

% quantum_interp - heuristic McCann interpolation-like
%
%   nu = quantum_interp(gamma, mu, m, d, options);
%
%   m is number of steps in [0,1] for the interpolation.
%
%   gamma is a matricial coupling
%   mu{1} and mu{2} are the marignals
%   nu{k} for k=1...m is the McCann-like interpolation
%
%   d should be either 1 (1-D case), 2 (2-D) or a generic cost matix (slow
%   mode).
%
%   options.sparse_mult controls the number of traveling diracs used.
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;

if size(d,1)>1
    c = d;
    d = 'slow';
elseif d~=1 && d~=2
    error('Only implemented in 1D and 2D (need to define the grid in higher dimensions)');
end

N = [size(mu{1},3) size(mu{2},3)];
if d==2
    n = floor(sqrt(N(1))); % width of the image
end

if isstruct(gamma)
    % already in sparse format
    I = gamma.i; J = gamma.j;
    gamma_ij = gamma.T;
else
    % convert in sparse format
    sparse_mult = getoptions(options, 'sparse_mult', 3);
    T = squeeze( trM(gamma) );
    v = sort(T(:), 'descend'); v = v(min(end,round(max(N)*sparse_mult))); % threshold
    [I,J] = ind2sub(N, find(T>v) );
    gamma_ij = gamma(:,:,find(T>v));
    mu0 = { sum(gamma,4) squeeze(sum(gamma,3)) };
end

% matched marginals
mu1 = { sparse_marginal(I,J,gamma_ij, N, 0), ...
        sparse_marginal(I,J,gamma_ij, N, 1) };
% mu1 = mu0;

% check
% plot([squeeze(trM(mu0{2})) squeeze(trM(mu1{2}))] );

% correction factor
A = {};
for k=1:2
    for i=1:N
        A{k}(:,:,i) = mu{k}(:,:,i)*pinv( mu1{k}(:,:,i) );
    end
end
% full interpolation
nu = {};

s = (0:N(1)-1)'/N(1); % 1D grid
if d==2
    s = (0:n-1)/n; % 2D grid
    [sY,sX] = meshgrid(s,s);
    s = [sX(:), sY(:)];
end


for k=1:m
    progressbar(k,m);
    t = (k-1)/(m-1);
    nu{k} = zeros(size(mu{1},1),size(mu{1},2),N(1),1);
    for i=1:length(I)
        % interpolated value
        G = ( (1-t)*A{1}(:,:,I(i)) + t*A{2}(:,:,J(i)) ) * ...
            gamma_ij(:,:,i);
        % interpolated position
        p = (1-t)*s(I(i),:) + t*s(J(i),:);
        % rounding on the grid
        switch d
            case 1
                %% 1D %%
                pi = floor(p*N)+1;
            case 2
                %% 2D %%
                pi = floor(p*n)+1;
                % turn it into a linear index
                pi = pi(1) + (pi(2)-1)*n;
            case 'slow'
                %% Generic, compute barycenter numerically %%
                [~,pi] = min( (1-t)*c(I(i),:) + t*c(J(i),:) );
        end
        % accumulate measure
        nu{k}(:,:,pi) = nu{k}(:,:,pi) + G;
    end
end

end
