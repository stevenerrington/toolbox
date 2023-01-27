function filenames = getDirFilenames(storageDir,label)

matFiles = dir(fullfile(storageDir,'*.mat'));
matFiles = {matFiles.name};
DMFCfiles = ~cellfun('isempty',strfind(matFiles,label));
filenames = matFiles(DMFCfiles);

end
