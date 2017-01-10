% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function img = align_imgs(path_source, path_target)

    [~,fn_t,~] = fileparts(path_target);
    [path,fn_s,~] = fileparts(path_source);
    out_fn = sprintf('%s\\%s_to_%s.png', path, fn_s, fn_t)
     
    source = double(imread(path_source))/255;
    target = double(imread(path_target))/255;

    [p_source,p_target] = cpselect(source, (max(0,target)),'Wait',true);

    tform = maketform('projective', p_source, p_target);

    img = imtransform(source, tform, 'bilinear', ...
        'size', size(target), ...
        'xdata',[0 size(source,2)], ...
        'ydata',[0 size(source,1)]);

    
    imwrite(img,out_fn);
    
%	imagec(img);

    
end