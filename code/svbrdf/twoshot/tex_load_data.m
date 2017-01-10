% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function [fdata, transport] = tex_load_data(datapath)
    [path, name] = standardize_path(datapath, 1)    
    load(strcat(path, 'out/fdata.mat'));
    fdata.path = path;
    fdata.name = name;
    
    transport = [];

    if nargout > 1 && exist(strcat(path, 'out/transport.mat'))

        load(strcat(path, 'out/transport.mat'));
        transport = res;
    end
end
