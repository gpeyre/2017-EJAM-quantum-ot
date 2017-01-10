% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

% Solve and output for all datasets

function SOLVE_ALL(datapath)
    
    if nargin == 0
        % If no path to data/ directory is supplied, use a default one
        % (this is convenient if you simply want to start your computations
        % from here by pressing F5)
        datapath = 'e:\texture_rel\test\data\';
    end

    % The default dataset list; feel free to delete this and use your own
    files = {
        'book_black',
        'book_bluegray',
        'book_brown',
        'book_leather_red',
        'book_pale',
        'book_red',
        'book_text',
        'book_wine',
        'book_worn',
        'book_yellow',
        'bread_crisp',
        'brick_red',
        'cardboard',
        'fabric_blue',
        'fabric_chair',
        'fabric_office',
        'fabric_orange',
        'fabric_table',
        'fabric_weave',
        'fabric_yellow',
        'fabric_zigzag',
        'laptop_bag',
        'leather_amber',
        'leather_antique',
        'leather_black',
        'leather_brown',
        'leather_dijon',
        'leather_white',
        'leather_wine',
        'metal_aniso',
        'metal_black',
        'metal_cyan',
        'metal_galvanized',
        'metal_gritty',
        'metal_rust',
        'metal_scratches',
        'paint_black',
        'paint_cyan',
        'paint_thick',
        'paint_white',
        'paper_crumpled',
        'plastic_blue',
        'plastic_cutting',
        'plastic_green',
        'plastic_red',
        'plastic_scratches',
        'plastic_weave',
        'plastic_weave_silver',
        'print',
        'rubber_parallel',
        'seed',
        'tape_silver',
        'tape_silver2',
        'tape_yellow',
        'tile_stair',
        'tinfoil',
        'tissue',
        'wallpaper2',
        'wallpaper3',
        'wallpaper4',
        'wallpaper5',
        'wallpaper6',
        'wallpaper_white',
        'wood_aniso',
        'wood_board',
        'wood_board_brown',
        'wood_dark',
        'wood_door',
        'wood_laminate',
        'wood_old',
        'wood_plain',
        'wood_shiny'};


	% if we want to skip some with a regular pattern (e.g. to run multiple 
    % Matlabs on the same set), then for example:
    % files = files(1:3:end)    % for Matlab #1
    % files = files(2:3:end)    % for Matlab #2
    % files = files(3:3:end)    % for Matlab #3
    % (note that running concurrent CUDA steps probably isn't a good idea
    % though)
    
    % Go through all the datasets and steps. It might be a good idea to
    % split this into two separate batches, first one doing all the CUDA
    % transports, and the second one doing the datafits in some concurrent
    % fashion.
    for fn = files'
        figure

        % Initial reflectance transport
        fdata = tex_load_data(char(strcat(datapath,fn)));     
        transport = tex_compute_transport_cuda(fdata);            
        
        % ... done, load the result and feed it to data fit step:
        [fdata,transport] = tex_load_data(char(strcat(datapath,fn)));     
        tex_alternate(fdata, transport, 'solution.mat');

        % ... and finally, reverse transport and output
        tex_reverse_output(strcat(datapath, char(fn)),'',true);
    end

end

