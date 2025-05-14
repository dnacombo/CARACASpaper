%%
setpath_ds003690
%% collect all files created by script_02_ds003690_precomp

load(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs')
clear ffs
for i_f = 1:numel(fs)
    fprintf('%s\n',myfileparts(fs(i_f).name,'f'));
    load(fullfile(dirout,fs(i_f).sub,sprintf('%s_%s_%s_FilesAndScoresList.mat',fs(i_f).sub, fs(i_f).task, fs(i_f).run)),'f')
    if not(strcmp(fs(i_f).name,f.name))
        error('mixing things up...')
    end
    ffs(i_f) = removefields(f,'MANUAL');
end
fs = ffs;
save(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs');
%% collect all ratings done by Pierre
load(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs');

Pierre = readtext('/network/iss/cenir/analyse/meeg/CARACAS/allfiles_ds3690_verif_Pierre.csv', ',');
Pierre = Pierre(:,1:5);
Paul = struct();
for i = 2:size(Pierre,1)
    for j = 1:size(Pierre,2)
        Paul(i-1).(Pierre{1,j}) = Pierre{i,j};
        if j == 5 && ischar(Paul(i-1).(Pierre{1,j}))
            Paul(i-1).(Pierre{1,j}) = str2double(strsplit(Paul(i-1).(Pierre{1,j}),';'));
        end
    end
end
Pierre = Paul;clear Paul

for i_f = 1:numel(fs)
    fprintf('%s\n',myfileparts(fs(i_f).name,'f'));
    P = flist_select(Pierre,'sub',['^' fs(i_f).sub '$'],'task',fs(i_f).task,'run',fs(i_f).run);
    if numel(P) ~= 1
        error('Ambiguous')
    end
    fs(i_f).MANUAL.rej = false(size(fs(i_f).CARACAS.rej));
    fs(i_f).MANUAL.rej(P.CardiComps) = 1;
end

save(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs');