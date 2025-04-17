setpath_ds003690


load(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs')
for i_f = 1:numel(fs)
    load(fullfile(dirout,fs(i_f).sub,sprintf('%s_FilesAndScoresList.mat',fs(i_f).sub)),'f')
    ffs(i_f) = f;
end
fs = ffs;
save(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs');
