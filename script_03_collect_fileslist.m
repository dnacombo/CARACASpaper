setpath_ds003690


load(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs')
clear ffs
for i_f = 1:numel(fs)
    fprintf('%s\n',myfileparts(fs(i_f).name,'f'));
    load(fullfile(dirout,fs(i_f).sub,sprintf('%s_%s_%s_FilesAndScoresList.mat',fs(i_f).sub, fs(i_f).task, fs(i_f).run)),'f')
    if not(strcmp(fs(i_f).name,f.name))
        error('mixing things up...')
    end
    ffs(i_f) = f;
end
fs = ffs;
save(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs');
