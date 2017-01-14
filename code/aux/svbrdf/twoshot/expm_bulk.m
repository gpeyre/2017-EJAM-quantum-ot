% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).


function E = expm_bulk(X)
    % see: Wikipedia talk page for matrix exponential
    k = size(X,2);
    Y = [X(1,:); X(3,:); X(3,:); X(2,:)];
    Y = reshape(Y, 2, 2, []);
    
    I = repmat([1 0; 0 1], [1 1 k]);
    TR = Y(1,1,:) + Y(2,2,:);
    
    T = Y;
    T(1,1,:) = T(1,1,:) - TR/2;
    T(2,2,:) = T(2,2,:) - TR/2;
    
    q = sqrt(- (T(1,1,:).*T(2,2,:) - T(2,1,:).*T(1,2,:)));
    seTR = (exp(TR/2));
    t1 = seTR .* cosh(q);
    t2 = seTR .* sinh(q) ./ q;
    
    E = repmat(t1, [2 2 1]).*I + repmat(t2, [2 2 1]) .* T;
end