
setpath_ds003690


load(fullfile(dirout,'AllFilesAndScoresList.mat'))

if ispc
    for i = 1:numel(fs)
        fields = {'name', 'compf','eegf'};
        for j = 1:numel(fields)
            fs(i).(fields{j}) = strrep(fs(i).(fields{j}),'/network','\');
            fs(i).(fields{j}) = strrep(fs(i).(fields{j}),'/','\');
        end
    end
end
%%

CORRrej = [];CARACASrej = [];SASICARACASrej = [];
for i = 1:numel(fs)
    CORRrej(i,:) = fs(i).CORR.rej;
    CARACASrej(i,:) = fs(i).CARACAS.rej;
    % CARACASrej(i,:) = [fs(i).CARACAS.meas.PQ] > 1/3;
    % SASICARACASrej(i,:) = fs(i).SASICARACAS.rej;
    MANUALrej(i,:) = fs(i).MANUAL.rej;
end
%%
figure(595);clf
set(gcf,'UserData',fs)
subplot(321)
hi = imagesc(MANUALrej');
title('MANUAL')
set(hi,"ButtonDownFcn",@plotit);
subplot(322);
hi = imagesc(CARACASrej');
title('CARACAS')
set(hi,"ButtonDownFcn",@plotit);

subplot(3,2,[3 6]);
toplot = NaN(size(MANUALrej));
toplot(MANUALrej & CARACASrej) = 1;
toplot(MANUALrej & ~ CARACASrej) = 2;
toplot(~MANUALrej & ~CARACASrej) = 3;
toplot(~MANUALrej & CARACASrej) = 4;
names = {'Hit','Miss','CR','FA'};

hi = imagesc(toplot');
vline(1:numel(fs),':k','pickableparts','none')
hline(1:size(toplot,2),':k','pickableparts','none')
set(gca,'Colormap',varycolor(4))
set(hi,"ButtonDownFcn",@plotit);
clim([0.5 4.5])
h = colorbar();
set(h,'YTick',1:4,'yticklabel',names)






function plotit(src,evt)

fs = get(gcf,'UserData');
p = get(gca,'CurrentPoint');
icmp = round(p(1,2));
ids = round(p(1,1));

load('layout.mat')

fig(21,4894);clf;
h = subplot(211);

comp = load(fs(ids).compf);
cfg = [];
cfg.layout = layout;
cfg.component = icmp;
cfg.channel = 'eeg';
cfg.figure = h;
ft_topoplotIC(cfg,comp)

title(['comp' num2str(icmp) ' ' fs(ids).sub,' ', fs(ids).task,' ', fs(ids).run], 'fontsize',16);


data = load(fs(ids).eegf);
EKGchan = chnb('ekg',data.label);
EKG = data.trial{1}(EKGchan,:);
EKG = EKG/range(EKG) * range(comp.trial{1}(icmp,:));

subplot(212)
plot(comp.time{1},comp.trial{1}(icmp,:))
hold on
plot(comp.time{1},EKG,'r:')

xl = xlim;yl = ylim;

toprint = removefields(fs(ids).CARACAS.meas(icmp),{'Ampl_var'});
fields = fieldnames(toprint);
threshs_min = [0 0 0 0 0 0 2 0 0 35];
threshs_max = [1/3 1/3 1/3 1/3 1/3 1/3 inf 1/3 1/3 90];
for i = 1:numel(fields)
    strtitle = [fields{i} ' = ' num2str(toprint.(fields{i}),2)];
    c = 'k';
    c = ifelse(toprint.(fields{i}) < threshs_min(i),'r',c);
    c = ifelse(toprint.(fields{i}) > threshs_max(i),'r',c);
    x = xl(1)+(diff(xl)/3)*rem(i-1,round(numel(fields)/3)) +1;
    y = yl(1)-(floor((i)/(numel(fields)/3))+1)*diff(yl)/20-diff(yl)/10;
    text(x,y,strtitle,'HorizontalAlignment','left','VerticalAlignment','baseline','FontSize',8, 'Interpreter','none', 'color',c)
end
title([ifelse(fs(ids).CARACAS.rej(icmp),'CARACAS ',''),ifelse(fs(ids).MANUAL.rej(icmp), ' MANUAL ' ,'')])

drawnow

cfgtmp = [];
if isfield(comp, 'sampleinfo')
    cfgtmp.trl = [1, comp.sampleinfo(end,2), 0];
else
    cfgtmp.trl = [1, sum(cellfun(@numel,comp.time)), 0]; % Start, End, Offset
end
comp = ft_redefinetrial(cfgtmp, comp);

comp.trial{1} = comp.trial{1}(icmp,:);
comp.label = comp.label(icmp);
comp.topolabel = comp.topolabel(icmp);

cfg = [];
cfg.method_chosen = 'absolute_threshold';   %'absolute_threshold' (method 1) or 'mean_std' (method 2)
cfg.plot_heart_IC = 0;                      % 1 or 0 (to plot the IC labelled as cardiac)
cfg.path_output = []; % path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
cfg.file_info = myfileparts(fs(ids).name,'f');% name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)
cfg.nb_IC_wanted = 3;                       % number of IC selected for each metric (kurtosis, skewness...) [default: 3, to select the top 3 IC for each metric]
cfg.bpm_min = 45;                           % expected heart beat per min, for sanity check [default: 45 and 90]
cfg.bpm_max = 90;
cfg.threshold_cond_IC_method1 = .5;         % minimum proportion of conditions that must be met in order that an IC could be considered as a potential heart IC [default: 0.5, so if method_chosen == 'absolute_threshold', an IC must be in the top 3 for at least 50% of the metrics]
cfg.threshold_std_method2 = 2.5;            % if method_chosen == 'mean_std', an IC will be considered as a potential heart IC if its proportion of conditions met (i.e., its score) is above mean(all_score) + threshold_std_method2 * std(all_score) [default: 2.5]
cfg.min_recording_duration_sec = 20;        % minimum duration (in sec) of the IC timecourse (default: 20]
cfg.mini_bouts_duration_for_SignalAmplRange = 10; % for sanity check (avoids false positive): the time course of a potential heart IC must be ~regular. The timecourse will be divided into mini-segments of this duration, and we will check that the amplitude between these mini-bouts is ~similar. [default: 10]
cfg.threshold_regularity_signal_minmax = 1.5; % For each mini-bout, the averaged signal amplitude is computed. The IC timecourse will be considered as irregular if: (max(Mean_Amp_minibout) - min(Mean_Amp_minibout)) / min(Mean_Amp_minibout) > threshold_regularity_signal_minmax [default: 1.5]

[isCardiac, meas] = CARACAS(cfg, comp);

ids

corr = fs(ids).CORR.c(icmp)

meas
isCardiac
% cfg = [];
% cfg.viewmode  = 'component';
% cfg.layout    = layout;
% ft_databrowser(cfg, comp);

end