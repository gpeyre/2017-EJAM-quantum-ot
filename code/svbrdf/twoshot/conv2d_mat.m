% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function A = conv2d_mat(K, s, border)
    h = s(1);
    w = s(2);

    % Create a flattened index array
    idx = reshape(1:(w*h), [h w]);
    
    % Expand it by the filter size, using supplied boundary conditions
    Kodd = mod(size(K)-1,2);
    Keven = size(K)-1-Kodd;
%    idx = padarray(idx, ceil([(size(K,1)-1)/2, (size(K,2)-1)/2]), border);
    idx = padarray(idx, Keven/2, border);
    idx = padarray(idx, Kodd, border, 'post');
    
    % Accumulate shifted and weighted neighbor according to K
    A = sparse(w*h, w*h);    
    for ki = 1:size(K,1)
        for kj = 1:size(K,2)
            M = spdiags(ones(w*h,1), 0, w*h, w*h);
            M = M(vec(idx(ki:h+ki-1, kj:w+kj-1)),:);
            A = A + K(ki,kj) * M;
        end
    end
end

