function nu = compute_quantum_interp(gamma, mu, m, d)

% compute_quantum_interp - heuristic McCann interpolation-like 
%
%   nu = compute_quantum_interp(gamma, mu, m, d);
%
%   gamma is a matricial coupling
%   mu{1} and mu{2} are the marignals
%   nu{k} for k=1...m is the McCann-like interpolation
%
%   d should be either 1 (1-D case), 2 (2-D) or a generic cost matix (slow
%   mode).
%
%   Copyright (c) 2016 Gabriel Peyre

if size(d,1)>1 
    c = d; 
    d = 'slow';
elseif d~=1 && d~=2
    error('Only implemented in 1D and 2D (need to define the grid in higher dimensions)');
end

N = size(gamma,3);
if d==2
    n = floor(sqrt(N)); % width of the image
end


% marginals
mu1 = { sum(gamma,4) squeeze(sum(gamma,3)) };
% correction factor
A = {};
for k=1:2
    for i=1:N
        A{k}(:,:,i) = mu{k}(:,:,i)*pinv( mu1{k}(:,:,i) );
    end
end
% full interpolation
nu = {};


s = (0:N-1)'/N; % 1D grid
if d==2
    s = (0:n-1)/n; % 2D grid
    [sY,sX] = meshgrid(s,s);
    s = [sX(:), sY(:)];
end

% threshold entries to speed up
a = trM(gamma); 
[I,J] = ind2sub(N, find(a>max(a(:))/100) );
for k=1:m
    t = (k-1)/(m-1);
    nu{k} = zeros(size(mu{1},1),size(mu{1},2),N,1);
    for i=1:length(I)
        % interpolated value
        G = ( (1-t)*A{1}(:,:,I(i)) + t*A{2}(:,:,J(i)) ) * ...
            gamma(:,:,I(i), J(i));
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