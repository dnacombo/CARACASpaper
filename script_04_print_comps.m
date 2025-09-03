function script_04_print_comps(i_f)

if not(exist('i_f','var'))
    i_f = 1;
end
which_data = '_SASICA';
setpath_ds003690

dont_plot_topos = 1;
trial2plot = 2;
decorate = 0;

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


for i_f = i_f%1:numel(fs)
    outpng = fullfile(dirout,'AllCompImages',[fs(i_f).sub,'_', fs(i_f).task,'_', fs(i_f).run, '_allComps_',ifelse(dont_plot_topos,'tc','topo'), '.png']);
    fprintf('%s\n',myfileparts(outpng,'f'));
    png = dir(outpng);
    compmat = dir(fs(i_f).compf);
    if 1%isempty(png) || png.datenum < compmat.datenum
        figure('position', [ 1         121        1920         988]);
        set(gcf,'WindowState','maximized')

        comp = load(fs(i_f).compf);
        data = load(fs(i_f).eegf);
        EKGchan = chnb('ekg',data.label);
        comp.cfg = [];
        %%
        for i_c = 1:numel(comp.label)
            if dont_plot_topos
                %%
                h = subplot(6,10,i_c);cla
                plot(comp.time{trial2plot},comp.trial{trial2plot}(i_c,:))
                hold on
                EKG = data.trial{trial2plot}(EKGchan,:);
                EKG = EKG/range(EKG) * range(comp.trial{trial2plot}(i_c,:));
                plot(comp.time{trial2plot},EKG,'r:')
                title(sprintf('Comp %d',i_c))
                xl = xlim;yl = ylim;
                if i_c == 1
                    text(xl(1),yl(2)+diff(yl)/2,[fs(i_f).sub,' ', fs(i_f).task,' ', fs(i_f).run], 'fontsize',18);
                end
                if decorate
                    toprint = removefields(fs(i_f).CARACAS.meas(i_c),{'Ampl_var'});
                    fields = fieldnames(toprint);
                    threshs = [];
                    threshs.PQ = [0 1/3];
                    threshs.QS = [0 1/3];
                    threshs.ST = [0 1/3];
                    threshs.PR = [0 1/3];
                    threshs.RT = [0 1/3];
                    threshs.PT = [0 1/3];
                    threshs.sk = [2 inf];
                    threshs.ku = [2 inf]; % Lower bound only used now, but keep range for visualization consistency
                    threshs.RR = [0 1/3];
                    threshs.Rampl = [0 1/3];
                    threshs.bpm = [35 90];
                    for i = 1:numel(fields)
                        strtitle = [fields{i} ' = ' num2str(toprint.(fields{i}),2)];
                        c = 'k';
                        c = ifelse(toprint.(fields{i}) < threshs.(fields{i})(1),'r',c);
                        c = ifelse(toprint.(fields{i}) > threshs.(fields{i})(2),'r',c);
                        text(xl(2),yl(2)-i*diff(yl)/10,strtitle,'HorizontalAlignment','left','VerticalAlignment','baseline','FontSize',6, 'Interpreter','none', 'color',c)
                    end
                    % text(xl(2),yl(2),strtitle,'HorizontalAlignment','left','VerticalAlignment','baseline','FontSize',6, 'Interpreter','none')
                    set(gca,'position',get(gca,'position') .* [1 1  .80 1])
                    if fs(i_f).CARACAS.rej(i_c)
                        set(gca,'XColor','r');
                    end
                    if isfield(fs(i_f), 'MANUAL') && fs(i_f).MANUAL.rej(i_c)
                        set(gca,'YColor','r');
                    end
                end
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
                plot(comp.time{trial2plot},comp.trial{trial2plot}(i_c,:))
            end
        end
        drawnow
        %%
        
        saveas(gcf,outpng)
        close(gcf)
    end
end