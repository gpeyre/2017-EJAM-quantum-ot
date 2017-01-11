dataDir = '/scratch/texture_data/';
textureTypes = strsplit(ls(data_dir));
nTextures = 5;

whichTextures = randperm(length(textureTypes),nTextures);

for i=1:nTextures
    type = textureTypes{whichTextures(i)}
    dir = [dataDir type '/'];
    [mat, S] = twoshot_load_dataset(dir);
    
    save(['svbrdf/data/' type '.mat'],'mat','S');
end