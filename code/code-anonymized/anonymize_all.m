%%
% Anonymize the files.

addpath('toolbox/');

str = input('Name to replace: ', 's');

rep_list = {...
    './', ...
    'toolbox/', ...
    'toolbox_anisotropic/', ...
    'toolbox_connections/', ...
    'toolbox_fast_marching/', ...
    'toolbox_geometry/', ...
    'toolbox_quantum/', ...
    };

for k=1:length(rep_list)
    anonymize(rep_list{k}, '.m', str);
end

 