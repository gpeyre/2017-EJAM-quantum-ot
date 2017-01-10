% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function Y = tex_filter_medioid_fast(X, n)
    nn = 2*n+1;

    S = size(X);
    
    Y = X;
    
    XP = padarray(X, [n n], 'both');
    
    D = repmat(XP, [1 1 1 nn nn]);
    for i = 1:nn
        for j = 1:nn
            D(:,:,:,i,j) = circshift(XP, [i-1-n j-1-n 0]);
        end
    end
    
    D = sum(sum(sqrt(sum((D - repmat(XP, [1 1 1 nn nn])).^2, 3)), 4), 5);
    
    D(1:n,:) = Inf;
    D(end-n+1:end,:) = Inf;
    D(:,1:n) = Inf;
    D(:,end-n+1:end) = Inf;
    
    
    
    for i = 1:S(1)
        for j = 1:S(2)
            w = D(i:i+nn-1, j:j+nn-1);  % (i,j) centered nn-window
            
            w = w(:);
            
            [~,m] = min(w);
            
            oi = mod(m-1,nn);
            oj = floor((m-1)/nn);
            
            Y(i,j,:) = X(max(1,min(S(1),i-n+oi)), max(1,min(S(2),j-n+oj)), :);
            
        end
    end
    
end