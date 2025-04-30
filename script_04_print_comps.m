setpath_ds003690

dont_plot_topos = 0;

mymkdir(fullfile(dirout,'AllCompImages'))

load(fullfile(dirout,'AllFilesAndScoresList.mat'))

load layout.mat
%%
s = keepfields(fs,{'sub','task','run'});
[s.id] = rep2struct(1:numel(s));
s = orderfields(s, {'id','sub','task','run'});
[s.CardiComps] = rep2struct('');
writetable(struct2table(s),'allfiles.csv');

%% 


for i_f = 1:numel(fs)
    outpng = fullfile(dirout,'AllCompImages',[fs(i_f).sub,'_', fs(i_f).task,'_', fs(i_f).run, '_allComps_',ifelse(dont_plot_topos,'tc','topo'), '.png']);
    fprintf('%s\n',myfileparts(outpng,'f'));
    png = dir(outpng);
    compmat = dir(fs(i_f).compf);
    if isempty(png) || png.datenum < compmat.datenum
        figure('position', [ 1         121        1920         988]);
        set(gcf,'WindowState','maximized')

        comp = load(fs(i_f).compf);
        comp.cfg = [];
        %%
        for i_c = 1:numel(comp.label)
            if dont_plot_topos

                h = subplot(6,10,i_c);
                plot(comp.time{1},comp.trial{1}(i_c,:))
                title(sprintf('Component %d',i_c))
            else
                h = subplot(12,10,i_c+10*(floor((i_c-1)/10)));
                cfg = [];
                cfg.layout = layout;
                cfg.component = i_c;
                cfg.channel = 'eeg';
                cfg.figure = h;
                cfg.comment = 'no';
                cfg.markersymbol = '.';
                ft_topoplotIC(cfg,comp)

                h = subplot(12,10,i_c+10*(floor((i_c-1)/10))+10);
                plot(comp.time{1},comp.trial{1}(i_c,:))
            end
        end
        drawnow
        saveas(gcf,outpng)
        close(gcf)
    end
end