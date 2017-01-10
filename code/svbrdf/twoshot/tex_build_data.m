% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function fdata = tex_build_data(datapath)

    % The tile resolution (M x M). However! The actual resolution will be
    % 2M x 2M for historical reasons. So M=96 will give the standard 
    % 192x192 tiles we used in the paper.
    M = 96;
    
    B = [12 12];

    [path, name] = standardize_path(datapath, 1)
    
    fdata = struct;
    
    fdata.path = path;
    fdata.name = name;

    fn = strcat(path,'guide.png');
    if ~exist(fn)
        fn = strcat(path,'guide.jpg');
        if ~exist(fn)
            fn = strcat(path,'guide.tiff');
        end
    end
    img = imread(fn);
    img_guide = (double(img)/255).^2.2;
    
    fn = strcat(path,'flash.png');
    if ~exist(fn)
        fn = strcat(path,'flash.jpg');
        if ~exist(fn)
            fn = strcat(path,'flash.tiff');
        end
    end
    img = imread(fn);
    img_flash = (double(img)/255).^2.2;
        
    S = size(img_guide);
    S = S(1:2);

    crop_c = [2 2];
    crop_s = min(S)*1-4;
    
    dispmode = 0;
    
    while true
        aspect = B(2)/B(1);
        wh = [crop_s crop_s*aspect];
        clf;
        if dispmode == 0
            imagec(img_guide)
        else
            imagec(img_flash)
        end
        
        axis equal;
        %axis xy;
        hold on;
        rectangle('Position',[crop_c(2) crop_c(1) wh(2) wh(1)],'EdgeColor','g','linewidth',2);
        
        stride = wh ./ B;
        
        for bi = 1:B(1)
            plot([crop_c(2) crop_c(2)+stride(1)*B(2)],[crop_c(1)+bi*stride(1) crop_c(1)+bi*stride(1)], '-k');
        end
        for bj = 1:B(2)
            plot([crop_c(2)+bj*stride(1) crop_c(2)+bj*stride(1)], [crop_c(1) crop_c(1)+stride(1)*B(1)], '-k');
        end
        
        [x,y,k] = ginput(1)
        if numel(k) == 0
            break;
        end
        
        switch k
            case 1
                crop_c = [y x];
            case 43
                crop_s = crop_s * 1.03;
            case 45
                crop_s = crop_s / 1.03;
            case 28
                B(2) = B(2) - 1;
            case 29
                B(2) = B(2) + 1;
            case 30
                B(1) = B(1) + 1;
            case 31
                B(1) = B(1) - 1;
            case 32
                dispmode = mod(dispmode+1,2);
        end
        % 28
           
        B
    end
    
    fdata.crop = struct;
    fdata.crop.s = crop_s;
    fdata.crop.c = crop_c;
    fdata.crop.s_img = size(img_guide);
    
    fdata.par.M = M;
    fdata.par.B = B;
    
    img_guide = img_guide(ceil(crop_c(1)):floor(crop_c(1)+crop_s), ...
                          ceil(crop_c(2)):floor(crop_c(2)+crop_s/B(1)*B(2)), :);
                      
    img_flash = img_flash(ceil(crop_c(1)):floor(crop_c(1)+crop_s), ...
                          ceil(crop_c(2)):floor(crop_c(2)+crop_s/B(1)*B(2)), :);

	S = size(img_guide);
    S = S(1:2);
    
    offs = S/2;
    rsize = S ./ B
    rsize = rsize(1);
    
    master_img = img_guide;
    fn = strcat(path,'master.png');
    if exist(fn)
        master_img = (double(imread(fn))/255).^2.2;
    end
	fn = strcat(path,'master.jpg');
    if exist(fn)
        master_img = (double(imread(fn))/255).^2.2;
    end
            
    while true
        disp('Click the master block position.');
        disp('Accept the current position by pressing Enter.');
        clf;
        imagec(master_img)
        axis equal;
        hold on;
        rectangle('Position', [offs(2) offs(1) rsize rsize],'EdgeColor','y');
        [x,y,k] = ginput(1)
        if numel(k) == 0
            break;
        end
        if k == 1
            offs = [y x];
        end
    end
    
    fdata.par.main_offs = offs ./ S;
    offs ./ S
    
    %{
    % Saving the originals is very wasteful, but we kept them just in case
    % (currently only the resolution of the original needs to be known, to
    % determine the crop factor in coordinate system construction)
    fdata.img_orig = struct;
    fdata.img_orig.guide = single(img_guide);   % Let's save some storage...
    fdata.img_orig.flash = single(img_flash);
    %}
    
    % Actually, let's use this:
    fdata.crop.img_orig_size = size(img_guide);
    
    I = B * M;
    
    fdata.img_std = struct;
    fdata.img_std.guide = imresize(img_guide,I,'bilinear');
    fdata.img_std.flash = imresize(img_flash,I,'bilinear');
    
    % img_double is what we're actually using (see the discussion on
    % resolution at the top)
    fdata.img_double = struct;
    fdata.img_double.guide = imresize(img_guide,I*2,'bilinear');
    fdata.img_double.flash = imresize(img_flash,I*2,'bilinear');
    
    clf;
    imagec(fdata.img_std.guide);
    axis equal;

    outdir = strcat(path, 'out/')
    if ~exist(outdir)
        mkdir(outdir)
    end
    
    
    
    save(strcat(outdir, 'fdata.mat'), 'fdata');
    
end
