% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function [img_diff, img_spec, img_normal, img_spec_shape, stats] = tex_encode_maps(l,g)

    S = size(l);

    % YUV conversion matrix
    YI = [1 0 1.13983; 1 -0.39465 -0.5806; 1 2.03211 0];
    % Specular color
    Scol =  max(0,YI* ([1;g(1);g(2)]));

    % Colored specular map, recall the exponent mapping (wasteful for storage...)
    img_spec = cat(3, Scol(1)*exp(l(:,:,6)), Scol(2)*exp(l(:,:,6)), Scol(3)*exp(l(:,:,6)));

    % NOTE: this code must correspond to the mappings in
    % tex_fit_reflectance and the model in tex_render

    % Get out of the logarithmic domain for diffuse albedo
    img_diff = exp(l(:,:,3:5));

    % Construct the SPD matrix corresponding to the glossiness/anisotropy
    % optimization variables
    M = reshape(l(:,:,7:9), [], 3)';
    M = expm_bulk(M);
    M(1,1,:) = M(1,1,:) + 0.001;
    M(2,2,:) = M(2,2,:) + 0.001;
    M = mat_inv_bulk2(M);
    
    M = permute(reshape(M, [4 S(1:2)]), [2 3 1]);

    %detM = (M(:,:,1).*M(:,:,4) - M(:,:,2).^2);

    % Now that we have the matrix, just store the upper triangle elements
    % as it's symmetric.
    img_spec_shape = M(:,:,[1 4 2]);
    
    % Normal map as unit vectors
    img_normal = cat(3, l(:,:,1:2), ones(size(l(:,:,1))));
    img_normal = img_normal ./ repmat(sqrt(sum(img_normal.^2,3)), [1 1 3]);

    return

    %% Now, some encoding to squeeze the values between [0,1] for typical file formats
    
    % Squeeze the X and Y axes of the normals (but not Z)
    img_normal(:,:,1:2) = 0.5 * img_normal(:,:,1:2) + 0.5;
    
    % Just assume we won't be having huge albedos, and there's enough
    % numerical accuracy.
    img_diff = 0.1*img_diff;
    img_spec = 0.1*img_spec;

    % Some bit more arbitrary mappings for the matrix coefficients,
    % designed to spread the data reasonably evenly within the range.
    img_spec_shape(:,:,1) = atan(img_spec_shape(:,:,1) / 100) / pi * 2; % [0,infty] -> [0,1]
    img_spec_shape(:,:,2) = atan(img_spec_shape(:,:,2) / 100) / pi * 2; % [0,infty] -> [0,1]
    img_spec_shape(:,:,3) = atan(img_spec_shape(:,:,3) / 100) / pi + 0.5; % [-infty,infty] -> [0,1]
    
    
    if nargout < 5
        return
    end
    
    ALL = [vec(img_spec_shape(:,:,1)) vec(img_spec_shape(:,:,2)) vec(img_spec_shape(:,:,3))];
    stats = [min(ALL); max(ALL); mean(ALL); std(ALL); prctile(ALL, 1); prctile(ALL,99); median(ALL);];
    
end