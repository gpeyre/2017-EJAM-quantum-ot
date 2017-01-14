% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function Z = wasserstein_project(X,Y,ntheta,iters)
    D = size(Y,1);

    % There is some capability for filtering data by having NaN's in it
    % (basically everything touched by a NaN is left out of the projection)
    % but this isn't used anymore.
	N0 = size(Y,2);
    temp = Y(:,~isnan(sum(Y,1)));
    if size(temp,2) == 0
        Z = X;
        return
    end
        
    Y = temp(:,count_resample(N0,size(temp,2)));

    for iter = 1:iters
        thetas = randn(D,ntheta);
        thetas = thetas ./ repmat(sqrt(sum(thetas.^2,1)), [D 1]);
       % [thetas,~] = qr(randn(D));
        
        H = zeros(D);
        for t = 1:ntheta
            H = H + thetas(:,t)*thetas(:,t)';
        end
        grads = H*X;
        
        for t = 1:ntheta
            theta = thetas(:,t);
            Xt = theta'*X;
            [~, Xto] = sort(Xt);
            Yt = theta'*Y;
            Yts = sort(Yt);
            Xdist = zeros(size(Yts));
            Xdist(Xto) = Yts;
            grads = grads - theta * Xdist;
            
        end
        
        X = X - pinv(H) * grads;
    end
    Z = X;
end