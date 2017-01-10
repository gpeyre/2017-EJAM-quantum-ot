% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function out = dice(img, M, mask)
    Simg = size(img);
    S = [M M Simg(3) Simg(1)/M Simg(2)/M];
    
    mask(mask < 0.001) = NaN;
    mask(~isnan(mask)) = 0;
    
    % Dice the flash image into blocks corresponding to 'recons'
    out = zeros(S);
    for bi = 1:S(4)
        for bj = 1:S(5)
            out(:,:,:,bi,bj) = (img((bi-1)*S(1)+1:bi*S(1), ...
                (bj-1)*S(2)+1:bj*S(2), :));
            
            if numel(mask) > 0
                out(:,:,:,bi,bj) = out(:,:,:,bi,bj) + repmat(mask(:,:,bi,bj), [1 1 3]);
            end
        end
    end

end