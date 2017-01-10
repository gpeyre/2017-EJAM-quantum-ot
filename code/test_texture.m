clear;

% put here the path where the data is
addpath('/Users/gpeyre/Dropbox/work/wasserstein/wasserstein-tensor-valued/matlab/svbrdf/data/');

addpath('svbrdf/');
addpath('svbrdf/twoshot/');
addpath('toolbox_quantum/');

rep = 'results/svbrdf/';
[~,~] = mkdir(rep);

%%
% Load a first material.

%load svbrdf/data/tile_stair.mat
load 'leather_white.mat';
whitemat = mat;

% render
viewing = twoshot_pointlight(S, [0.3 0.5 1], 1);
img = twoshot_render(mat, viewing);
imagec(img); title('leather\_white');
saveas(gcf,[rep 'leather_white.png']);

%%
% Load another material (replace with whatever path you are using)

load 'metal_black.mat';

[e1,e2,l1,l2] = tensor_eigendecomp(permute(mat.gloss, [3 4 1 2]));
A = (l1-l2)./(l1+l2); % anisotropy
E = l1+l2; % energy

% Render
viewing = twoshot_pointlight(S, [0.3 0.5 1], 1);
img = twoshot_render(mat, viewing);
figure; imagec(img); title('metal\_black');
saveas(gcf,[rep 'metal_black.png']);

%%
% Now, play with the gloss

mat2 = mat;
mat2.gloss = whitemat.gloss(1:S(1),1:S(2),:,:);
img = twoshot_render(mat2, viewing);
figure; imagec(img); title('metal\_black except leather\_white gloss');
saveas(gcf, [rep 'metal_black_with_leather_white_gloss.png']);

mat2 = mat;
mat2.gloss = .5*whitemat.gloss(1:S(1),1:S(2),:,:) + .5*mat.gloss;
img = twoshot_render(mat2, viewing);
figure; imagec(img); title('Average gloss');
saveas(gcf, [rep 'metal_black_with_average_gloss.png']);