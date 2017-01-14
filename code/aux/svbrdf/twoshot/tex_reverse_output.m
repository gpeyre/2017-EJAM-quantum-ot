% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function tex_reverse_output(datapath, ver, render)

    [path, name] = standardize_path(datapath, 1)
    load(sprintf('%s/out/solution%s.mat', path, ver));
    
    [fdata,transport] = tex_load_data(path);
    outdir = sprintf('%s/out/reverse%s', path, ver)
    try
        mkdir(outdir);
    catch err
        err
    end

    %{
    % Statistics collection, we used this to determine the file encoding...
    global stats_all;
    [~,~,~,~,stats] = tex_encode_maps(reshape(solution.step{2}.vars.l, [192 192 9]), solution.step{2}.vars.g); 
    stats_all = cat(3, stats_all, stats);
    size(stats_all)
    return
    %}
    
    %{
    % Just encode and show the master tile
    [img_diff, img_spec, img_normal, img_spec_shape] = tex_encode_maps(reshape(solution.step{2}.vars.l, [192 192 9]), solution.step{2}.vars.g); 
    clf
    subplot(2,3,1);
    imagec(img_diff);
    subplot(2,3,2);
    imagec(img_spec*0.1);
    subplot(2,3,3);
    imagec(img_normal*0.5+0.5);
    subplot(2,3,4);
    imagec(img_spec_shape);
    drawnow;
    return
    %}
    
    if 1
        R = tex_reverse_transport(fdata, transport, solution);
        save(strcat(outdir,'/reverse.mat'), 'R');
    else
        load(strcat(outdir,'/reverse.mat'));
    end
    
    S = size(R);

    STEP = 2;
    l = R;
    g = solution.step{STEP}.vars.g; 

    % Center the normal map (not a big deal but why not)
    l(:,:,1) = l(:,:,1) - mean(vec(l(:,:,1)));
    l(:,:,2) = l(:,:,2) - mean(vec(l(:,:,2)));

    [img_diff, img_spec, img_normal, img_spec_shape] = tex_encode_maps(l,g);
    
    flip = @(x) flipdim(x,1);

    
    write_pfm(flip(img_diff), strcat(outdir, '/map_diff.pfm'), 0);
    write_pfm(flip(img_spec), strcat(outdir, '/map_spec.pfm'), 0);
    write_pfm(flip(img_normal), strcat(outdir, '/map_normal.pfm'), 0);
    write_pfm(flip(img_spec_shape), strcat(outdir, '/map_spec_shape.pfm'), 0);
    

    clf
    subplot(2,3,1);
    imagec(img_diff);
    subplot(2,3,2);
    imagec(img_spec);
    subplot(2,3,3);
    imagec(img_normal.^2);
    subplot(2,3,4);
    imagec(img_spec_shape);

    drawnow;
    
    % Global variables
    fname = strcat(outdir, '/map_params.dat');
    fid = fopen(fname, 'wb');
    fprintf(fid, sprintf('%f %f', exp(g(4))+0.5, exp(g(3))+0.3));
    fclose(fid);    

    % Make a little rendered video if requested
    if render
        tex_render_rotate(R, solution.step{STEP}.vars.g, 0.5, strcat(outdir, '/rotate.avi'));
    end
end