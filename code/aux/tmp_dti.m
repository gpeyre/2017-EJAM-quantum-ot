
%% DATASETS (to download and put somewhere on your drive) %%
% - Example_data_for_ExploreDTI.mat from 
%     http://www.exploredti.com/exampledataset.htm
% - gk2-rcc-mask.raw from 
%     http://www.sci.utah.edu/~gk/DTI-data/
% - http://www.humanconnectome.org/data/
%       db.humanconnectome.org/data/projects/MGH_DIFF
% - http://www.mri-resource.kennedykrieger.org/databases
% - http://www.nitrc.org/projects/multimodal

% - http://sciviscontest.ieeevis.org/2010/


addpath('toolbox/');
addpath('toolbox_quantum/');
% change this to the path on your hard drive
addpath('/Users/gpeyre/Documents/data-dti/');

load('Example_data_for_ExploreDTI.mat');

trM = @(T)squeeze( T(1,1,:,:) + T(2,2,:,:) + T(3,3,:,:) );

% chosen depth
z = 30;

N = size( DT{1} ); N = N(1:2);

% ROI
selx = 14:111;
sely = 26:103;
N = [length(selx), length(sely)];

T = [];
T(:,:,1,1) = DT{1}(selx,sely,z);
T(:,:,1,2) = DT{2}(selx,sely,z);
T(:,:,1,3) = DT{3}(selx,sely,z);
T(:,:,2,2) = DT{4}(selx,sely,z);
T(:,:,2,3) = DT{5}(selx,sely,z);
T(:,:,3,3) = DT{6}(selx,sely,z);
% symmetrize
T(:,:,1,2) = T(:,:,1,2);
T(:,:,1,3) = T(:,:,3,1);
T(:,:,2,3) = T(:,:,3,2);
T = permute(T, [3 4 1 2]);

[Va,Ve] = Eigens3x3( reshape(T, [3 3 prod(N)]) );
Va = reshape(Va, [3 N]); Va = permute(Va, [2 3 1]);


% display trace, should be >=0
normalize = @(x)x/max(x(:));

anis = @(u,v)(u-v)./(u+v+1e-10);
Func = {@(x)sum(x,3), @(x)min(x,[],3), @(x)max(x,[],3), ...
       @(x)anis(max(x,[],3),min(x,[],3)) };
lgd = {'sum', 'min', 'max', 'aniso'};
clf;
for i=1:length(Func)
subplot(2,2,i);
imagesc(normalize(Func{i}(Va)));
axis image; axis off; colorbar;
colormap parula(256);
title(lgd{i});
end
