% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).


function res = brief_blur(img, WW, N, B)
    H = size(img,1);
    W = size(img,2);
    S = 2*WW+1;

    if B > 0
        G = fspecial('gaussian',ceil([B B]*4),B);
        img = imfilter(img,G,'same');    
    end
    
    X = randn(2,32*N) * 1/5 * S;
    Y = randn(2,32*N) * 1/5 * S;
    
    X = round(max(-WW,min(WW,X)));
    Y = round(max(-WW,min(WW,Y)));
%{    
    hold on;
    
    for i = 1:32*N
        plot(X(:,i),Y(:,i),'-');
    end
  %}  
    
    res = uint32(zeros(H,W,N));
    
    for i = 1:N
        for j = 1:32
            t = (i-1)*32+j;
            
            res(:,:,i) = res(:,:,i) + uint32(bitshift(1,j-1)) * ...
                uint32(circshift(img,[X(1,t) X(2,t)]) < circshift(img,[Y(1,t) Y(2,t)]));
            
          %  imagesc(uint32(circshift(img,[X(1,t) X(2,t)]) < circshift(img,[Y(1,t) Y(2,t)])));
           % drawnow;
        end
    end
end