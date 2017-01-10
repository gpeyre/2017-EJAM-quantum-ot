% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function out = heeger_bergen_multi_new(in, initial, iters, varargin)
    D = size(in,3);
    W = size(in,1);
    
    N_lev = log(W)/log(2)+2;
    
    if nargin > 3
    %    max_lev = varargin{1};
         N_lev = varargin{1};
    end
    
    witers = 5;
    wthetas = 5;
    
    pyr_in = buildSpyr2(in, N_lev, 'sp5.mat');
    
    N_or = numel(pyr_in{2});

    out = initial;
    flat = @(x) x(:);
    
    for iter = 1:iters
        pyr_out = buildSpyr2(out,N_lev,'sp5.mat');
        pyr_mod = pyr_out;
        

        pyr_mod{1}(2:end-1,2:end-1,:) = wasserstein_project_2d(pyr_out{1}(2:end-1,2:end-1,:),pyr_in{1}(2:end-1,2:end-1,:),wthetas,witers);
        for l = 2:N_lev-1
            for o = 1:N_or
                pyr_mod{l}{o}(2:end-1,2:end-1,:) = wasserstein_project_2d(pyr_out{l}{o}(2:end-1,2:end-1,:),pyr_in{l}{o}(2:end-1,2:end-1,:),wthetas,witers);
            end
        end
           
        if nargin <= 3
        	pyr_mod{end} = wasserstein_project_2d(pyr_out{end},pyr_in{end},wthetas,witers);
        end
            
        out = reconSpyr2(pyr_mod,'sp5.mat');

        out = wasserstein_project_2d(out,in,wthetas,witers);
    end
end


function Z = wasserstein_project_2d(X,Y,ntheta,iters)
    D = size(X,3);
    Z = wasserstein_project(reshape(permute(X,[D 1:(D-1)]),D,[]), ...
            reshape(permute(Y,[D 1:(D-1)]),D,[]),ntheta,iters);
	Z = reshape(Z', size(X));
end

