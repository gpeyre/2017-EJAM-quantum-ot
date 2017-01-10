clear

%load svbrdf/data/tile_stair.mat
load svbrdf/data/leather_white.mat
viewing = twoshot_pointlight(S, [0.3 0.5 1], 1);

whitemat = mat;

img = twoshot_render(mat, viewing);
imagec(img); title('leather\_white');
saveas(gcf,'leather_white.png');

% Load a material (replace with whatever path you are using)
load svbrdf/data/metal_black.mat
viewing = twoshot_pointlight(S, [0.3 0.5 1], 1);

% Render
img = twoshot_render(mat, viewing);
figure; imagec(img); title('metal\_black');
saveas(gcf,'metal_black.png');

% Now, play with the gloss
mat2 = mat;
mat2.gloss = whitemat.gloss(1:S(1),1:S(2),:,:);
img = twoshot_render(mat2, viewing);
figure; imagec(img); title('metal\_black except leather\_white gloss');
saveas(gcf,'metal_black_with_leather_white_gloss.png');

mat2 = mat;
mat2.gloss = .5*whitemat.gloss(1:S(1),1:S(2),:,:) + .5*mat.gloss;
img = twoshot_render(mat2, viewing);
figure; imagec(img); title('Average gloss');
saveas(gcf,'metal_black_with_average_gloss.png');