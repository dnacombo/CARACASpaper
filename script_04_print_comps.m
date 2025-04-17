setpath_ds003690


load(fullfile(dirout,'AllFilesAndScoresList.mat'))

load layout.mat

%% 
figure(4048);
set(gcf,'WindowState','maximized')

for i_f = 1:numel(fs)
    clf
    comp = load(fs(i_f).compf);
    %%
    for i_c = 1:numel(comp.label)

        h = subplot(12,10,i_c+10*(floor((i_c-1)/10)));
        cfg = [];
        cfg.layout = layout;
        cfg.component = i_c;
        cfg.channel = 'eeg';
        cfg.figure = h;
        cfg.comment = 'no';
        cfg.checkstring = 'no';
        ft_topoplotIC(cfg,comp)

        h = subplot(12,10,i_c+10*(floor((i_c-1)/10))+10);
        plot(comp.time{1},comp.trial{1}(i_c,:))
    end
    drawnow
    saveas(gcf,fullfile(dirout,fs(i_f).sub, fs(i_f).mod,[fs(i_f).sub,'_', fs(i_f).task,'_', fs(i_f).run, '_allComps.png']))
end