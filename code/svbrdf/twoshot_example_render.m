function twoshot_example_render(datapath)
    % Load a material (replace with whatever path you are using)
    [mat, S] = twoshot_load_dataset(datapath);

    % Build a basic coordinate system and a point light illumination
    % condition
	viewing = twoshot_pointlight(S, [-0.3 0.5 1], 1);
    
    % Render
    img = twoshot_render(mat, viewing);
    
    imagec(img);
end
