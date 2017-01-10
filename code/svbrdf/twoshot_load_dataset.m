% Loads a dataset from the same output pfm files that the C++ renderer uses

function [mat, S] = twoshot_load_dataset(path_in)
    path_in = standardize_path(path_in, 1);
    
    img_diff = read_pfm(strcat(path_in, 'out/reverse/map_diff.pfm'), true);
    img_spec = read_pfm(strcat(path_in, 'out/reverse/map_spec.pfm'), true);
    img_normal = read_pfm(strcat(path_in, 'out/reverse/map_normal.pfm'), true);
    img_spec_shape = read_pfm(strcat(path_in, 'out/reverse/map_spec_shape.pfm'), true);

    S = size(img_diff);
    S = S(1:2);
    
    M = cat(3, img_spec_shape(:,:,1), ...
               img_spec_shape(:,:,3), ...
               img_spec_shape(:,:,3), ...
               img_spec_shape(:,:,2));
	M = reshape(M, [S 2 2]);
        
    mat = struct;
    mat.alb_diff = img_diff;
    mat.alb_spec = img_spec;
    mat.normal = img_normal;
    mat.gloss = M;
    
    fid = fopen(strcat(path_in, 'out/reverse/map_params.dat') ,'r');
    B=textscan(fid,'%f ', 'delimiter',',');
    C=cell2mat(B);
    fclose(fid);
    mat.alpha = C(1);
end

