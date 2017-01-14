% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function coords = tex_form_coordinates(fdata, hfov, res_full)
    % hfov is the width of the image in a coordinate system where the
    % camera is perpendicularly viewing the surface from distance of 1. For
    % iPhone 5 it seems to be about 1.09 but varies based on focus
    % distance and such.
    %
    % res_full is the original size of the uncropped image (2448 x 3264 for
    % iPhone 5)

    % This is all very clumsy but whatever. There's also some messy legacy
    % stuff wrt. resolution (hence lots of doubling and halving...)
        
    M = fdata.par.M;
    
    res_crop = [];
    
    if isfield(fdata.crop, 'img_orig_size')
        res_crop = fdata.crop.img_orig_size;
    else
        res_crop = size(fdata.img_orig.flash);
    end
	res_crop = res_crop(1:2)
    
    res_std = size(fdata.img_double.flash) / 2;
    res_std = res_std(1:2)
    
    crop = res_crop ./ res_full;
    
    asp_full = res_full(2)/res_full(1);
    asp_crop = res_crop(2)/res_crop(1);
    cJ = repmat(1:res_std(2), [res_std(1) 1]);
    cI = repmat((1:res_std(1))', [1 res_std(2)]);
    
    [center_i,center_j] = img_max2(fdata.img_double.flash)
    center_i = (center_i)/2;
    center_j = (center_j)/2;
%    plot(j,i,'ro');
    
    % center the coordinate system at the highlight
    cI = -(cI - center_i);  % flip y/i-axis to get right-handed coordinates
    cJ = cJ - center_j;
    
    cI = cI / res_std(1) * crop(1) / asp_full * hfov;
    cJ = cJ / res_std(2) * crop(2) * hfov;
    
    cI = cI(M/2:M:end, M/2:M:end);
    cJ = cJ(M/2:M:end, M/2:M:end);

    
    fdata.par.B
      
    coords = cat(3, cJ, cI);    % flip, to get X,Y
end


function [i,j] = img_max2(img)
    % Blur the image slightly, find the centroid of the 5% of the brightest
    % pixels.
    H = size(img,1);
    W = size(img,2);
    I = repmat((1:H)', [1 W]);
    J = repmat((1:W), [H 1]);

    img = sum(img,3);
    

    GG = exp(-linspace(-2,2,ceil(min(size(img))/32)).^2);
    GG = GG / sum(GG);

    img2 = imfilter(img,GG,'same','symmetric');    
    img2 = imfilter(img2,GG','same','symmetric');    
 
    img2 = img2 - min(img2(:));

    img2 = img2 + 0.001*randn(size(img2));  % just to randomly scramble ordering in saturated pixels

    ss = sort(img2(:));
    cut = ss(round(numel(img2)*0.95));
    img_thresh = (img2 > cut);
 
    s = sum(sum(img_thresh));
    i = sum(sum(img_thresh .* I)) / s;
    j = sum(sum(img_thresh .* J)) / s;

end
