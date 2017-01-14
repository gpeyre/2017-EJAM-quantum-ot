% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function A = sparse_clip_boundary(s, b)
    h = s(1);
    w = s(2);
    
    idx = reshape(1:(h*w), [h w]);
    idx = idx(1:end-b(1), 1:end-b(2));
    
    A = spdiags(ones(h*w,1), 0, h*w, h*w);
    A = A(idx(:),:); 
end
