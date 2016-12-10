function u = lse(K)

% lse - Log-sum-exp
%
%   u = lse(K);
%
% compute u_i=log o sum_j o exp(K_ij)
%
%   You can set the global variable lse_stab. 
%   Use lse_stab==1 to use stabilized version using the shifting trick.
%   This should not be needed here, because the Sinkhorn iterates are
%   written in a stabilized way.
%
%   Copyright (c) 2016 Gabriel Peyre

global lse_stab;
if isempty(lse_stab)
    lse_stab = 0;
end

switch lse_stab
    case 0
        u = logM( sum(expM(K),4) );
    case 1
        d = size(K,1); % dimension
        if d~=2
            warning('only implemented in 2D');
        end
        t = ( K(1,1,:,:)+K(2,2,:,:) )/d;
        t = max(t,[],4);
        a = tensor_diag(t,t);
        A = repmat(a, [1 1 1 size(K,4)]);
        u = a+logM(sum(expM(K-A),4) );
end

end