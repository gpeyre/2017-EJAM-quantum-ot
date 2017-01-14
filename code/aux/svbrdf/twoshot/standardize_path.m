% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

% Just a quick and dirty function for sandardizing \'s and /'s and 
% whatever in file paths

function [path, name] = standardize_path(path, check)

    path = regexprep(path, '\\', '/');  % \'s into /'s
    path = regexprep(path, '/$', '');       % remove trailing /
    [~,name,~] = fileparts(path);
    path = strcat(path, '/');               % ... and reintroduce
    
    if check && exist(path) ~= 7
        error(sprintf('Directory %s not found', path'));
    end
    
end