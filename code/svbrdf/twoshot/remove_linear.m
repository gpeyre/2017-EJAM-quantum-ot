% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function img = remove_linear(img)
    s = [size(img) 1];
    
    g1 = repmat((1:s(1))', [1 s(2)]);
    g1 = g1 - mean(g1(:));
    g1 = g1 / norm(g1, 'fro');
    %g1 = g1(:);
    
    g2 = repmat((1:s(2)), [s(1) 1]);
    g2 = g2 - mean(g2(:));
    g2 = g2 / norm(g2, 'fro');
    %g2 = g2(:);
    
    g12 = g1.*g2;
    g12 = g12 / norm(g12);
    g11 = g1.^2;
    g11 = g11 - mean(g11(:));
    g11 = g11 / norm(g11);
    g22 = g2.^2;
    g22 = g22 - mean(g22(:));
    g22 = g22 / norm(g22);
    

    img = reshape(img, s(1), s(2), prod(s(3:end)));
    
    for c = 1:size(img,3)
        img(:,:,c) = img(:,:,c) - sum(vec(img(:,:,c).*g1)) * g1;
        img(:,:,c) = img(:,:,c) - sum(vec(img(:,:,c).*g2)) * g2;
        %{
        img(:,:,c) = img(:,:,c) - sum(vec(img(:,:,c).*g11)) * g11;
        img(:,:,c) = img(:,:,c) - sum(vec(img(:,:,c).*g12)) * g12;
        img(:,:,c) = img(:,:,c) - sum(vec(img(:,:,c).*g22)) * g22;
        %}
    end
    
    img = reshape(img, s);
    
    
end

