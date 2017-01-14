% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function I = imagen(V, a)    
    N = cat(3, V, ones(size(V,1), size(V,2)));
    I = N ./ repmat(sqrt(sum(N.^2, 3)), [1 1 3]);
    
    I(:,:,1:2) = I(:,:,1:2)+0.5;
    
    imagesc(max(0,min(1,I)));
end