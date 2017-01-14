% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function PP = tex_filter_transport(P, fs)

    P = double(P);
    
    S = size(P);

    I = repmat((1:S(1))', [1 S(2)]);
    J = repmat(1:S(2), [S(1) 1]);
    
    SI = I(P)-I;
    SJ = J(P)-J;
    
    % We are now in "offset" domain where "rigid" transport is seen as
    % constant values

    % Now filter to get flatter regions
    SISJ = tex_filter_medioid_fast(cat(3,SI,SJ), fs);
    SI = SISJ(:,:,1);    
    SJ = SISJ(:,:,2);    
    

    % Let's go back to index domain
    
    IP = SI+I;
    JP = SJ+J;

    IP = max(1, min(S(1), IP)); % ... can this even happen?
    JP = max(1, min(S(2), JP));
    
    
    PP = S(1)*(JP-1) + IP;
    
    subplot(3,2,5);
    imagesc(P);
    subplot(3,2,6);
    imagesc(PP);

    norm(P-PP,'fro');
    
end