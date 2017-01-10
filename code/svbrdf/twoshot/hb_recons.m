% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function out = hb_recons(source, target, Nlev)

    if nargin < 3
        Nlev = 5;
    end
        

    S = size(source);
    
    out = zeros(S);
    
    clf;
    
    % Loop through all the tiles and do the H-B
    for bi = 1:S(4)
        for bj = 1:S(5)
            [bi bj]
            out(:,:,:,bi,bj) = heeger_bergen_multi_new(target(:,:,:,bi,bj), source(:,:,:,bi,bj), 4, Nlev);
            
            subplot(2,2,1);
            imagec(out(:,:,:,bi,bj));
            subplot(2,2,2);
            imagec(target(:,:,:,bi,bj));
            subplot(2,2,3);
            imagec(source(:,:,:,bi,bj));
            drawnow;
        end
    end

end