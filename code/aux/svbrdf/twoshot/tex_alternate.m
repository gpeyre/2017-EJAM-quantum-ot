% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function tex_alternate(fdata, transport, solname)

    solution = struct;

    % Current estimate of the lumitexels (from gradient-domain initial 
    % transport)
    recons = transport.recons1;
    S = size(recons);

    % Solution utput directory
    fnbase = strcat(fdata.path, 'out/')
    
    % XXX! Hardcoded iphone 5 resolution and FOV! These should rather be
    % introduced somewhere in the data creation.
    coords = tex_form_coordinates(fdata, 1.09, [2448 3264]);
    

	[vars, rend] = tex_fit_reflectance(recons, fdata, coords, [], 7);
    
    solution.step{1}.vars = vars;
    solution.step{1}.rend = rend;
    
    save(strcat(fnbase, solname), 'solution');

    
    %load(strcat(fnbase, solname));
    
    
    % Add a tiny bit of identical noise to all tiles, just to encourage 
    % some variety
    noise = randn(S(1:2));
    %{
    noise = spatialPattern(S(1:2), -1);
    noise = noise - mean(noise(:));
    noise = noise / std(noise(:));
    %}
    rend_pre_hb = solution.step{1}.rend; % + 0.03*std(solution.step{1}.rend(:))  * repmat(noise, [1 1 S(3:5)]);    
    for i = 1:S(4)
        for j = 1:S(5)
            for c = 1:S(3)
               rend_pre_hb(:,:,c,i,j) = max(0,rend_pre_hb(:,:,c,i,j) + 0.03*std(vec(rend_pre_hb(:,:,c,i,j)))*noise);
            end
        end
    end

    % Assume any "multiplicative linear gradient" on the tile is caused by
    % the fact that the light direction may not be entirely constant within
    % the tile; remove it.
    rend_hb = hb_recons(rend_pre_hb, remove_linear_log( dice(fdata.img_double.flash,S(1),[]) ));
    
    solution.step{1}.rend_hb = rend_hb;
    
    % Second fitting, this time using the previous result as the initial
    % guess
    [vars2, rend2] = tex_fit_reflectance(rend_hb, fdata, coords, vars, 5);

    solution.step{2}.vars = vars2;
    solution.step{2}.rend = rend2;
    
    save(strcat(fnbase, solname), 'solution');

end