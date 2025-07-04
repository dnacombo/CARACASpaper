

% clear all
% close all

%% !!!!! USER CHOOSE PARAMETERS !!!!!

which_data = '_SASICA'; %''; %'_30comp';%   //

nb_IC_of_interest = 'all'; % 'all' // '10';%

which_col_in_csv = '_sure'; % '' // '_sure' // 'rej_noisy'

compare_truth = 'MANUAL';
compare_with = 'CARACAS'; %'SASICARACAS'; % 'CARACAS'; % 'CORR';


%% Path and load data

setpath_ds003690

load(fullfile(dirout,'AllFilesAndScoresList.mat'))

if ispc
    for i_f = 1:numel(fs)
        fields = {'name', 'compf','eegf'};
        for j = 1:numel(fields)
            fs(i_f).(fields{j}) = strrep(fs(i_f).(fields{j}),'/network','\');
            fs(i_f).(fields{j}) = strrep(fs(i_f).(fields{j}),'/','\');
        end
    end
end


% %% Plot hist to find thresholds
%
%
% all_metric_of_interest = [];
% for i_f = 1:numel(fs)
%     all_metric_of_interest(i_f,:) = [fs(i_f).CARACAS.meas.bpm]; % .sk // RR // Rampl // bpm // ku
%     allrej(i_f,:)= [fs(i_f).MANUAL.rej];
% end
%
% nu_metric = all_metric_of_interest(allrej);
% nunu_metric = all_metric_of_interest(~allrej);
% figure;
% h1 = histogram(nu_metric, 'Normalization', 'probability');
% hold on
% h2 = histogram(nunu_metric, 'Normalization', 'probability'); % 'NumBins', 1000,
% hold on;
% xline(35, 'r--', 'LineWidth', 2);
% xline(95, 'r--', 'LineWidth', 2);
% xlabel('bpm')
% ylabel('% of cardiac and non-cardiac IC')
% legend([h1, h2], {'Cardiac IC', 'Non-cardiac IC'}, 'Location', 'best');
% set(gca, 'FontSize', 38);
%


%% Overwrite CARACAS to test specific thresholds for each variable


%%% Adapt thresholds here
cfg_SASICA = SASICA('getdefs');
cfg_CARACAS = cfg_SASICA.CARACAS;
cfg_CARACAS.enable = 1;

truthrej = NaN(numel(fs),size(fs(1).CORR.rej,2));
withrej = truthrej;

for i_f = 1:numel(fs)
    if strcmp(compare_truth,'MANUAL')
        truthrej(i_f,:) = fs(i_f).(compare_truth).(['rej' which_col_in_csv]);
    elseif ~isempty(regexp(compare_truth,'CARACAS', 'once'))
        truthrej(i_f,:) = CARACAS_rethresh(truthrej(i_f,:),fs(i_f).(compare_truth),cfg_CARACAS);
    else
        truthrej(i_f,:) = fs(i_f).(compare_truth).rej;
    end


    if ~isempty(regexp(compare_with,'CARACAS', 'once'))

        withrej(i_f,:) = CARACAS_rethresh(withrej(i_f,:),fs(i_f).(compare_with),cfg_CARACAS);
    else
        withrej(i_f,:) = fs(i_f).(compare_with).rej;
    end
end



%% Plot and performances

% Plot the two figures on the top (MANUAL and CARACAS)
figure(595);clf
set(gcf,'UserData',struct('fs',fs, 'compare_with', compare_with, 'compare_truth',compare_truth))
subplot(321)
hi = imagesc(truthrej');
title(compare_truth)
set(hi,"ButtonDownFcn",@plotit);
subplot(322);
hi = imagesc(withrej');
title(compare_with)
set(hi,"ButtonDownFcn",@plotit);

subplot(3,2,[3 6]);
toplot = NaN(size(truthrej));
toplot(truthrej & withrej) = 1;     % Hit
toplot(truthrej & ~ withrej) = 2;   % Miss
toplot(~truthrej & ~withrej) = 3;   % CR
toplot(~truthrej & withrej) = 4;    % FA
names = {'Hit','Miss','CR','FA'};


% Select only a subset of IC
if ~ strcmp(nb_IC_of_interest, 'all')
    toplot = toplot(:,1:str2num(nb_IC_of_interest));
end


% Plot overlaping the two (with Hit, Miss, FA, CR)
hi = imagesc(toplot');
vline(1:numel(fs),':k','pickableparts','none')
hline(1:size(toplot,2),':k','pickableparts','none')
set(gca,'Colormap',varycolor(4))
set(hi,"ButtonDownFcn",@plotit);
clim([0.5 4.5])
h = colorbar();
set(h,'YTick',1:4,'yticklabel',names)


%%%% Compute perf metrics
% Sensitivity = Hit / (Hit + Miss)
Sensitivity = sum(toplot == 1) / (sum(toplot== 1) + sum(toplot == 2))
% Specificity = CR  / (CR + FA)
Specificity = sum(toplot == 3) / (sum(toplot== 3) + sum(toplot == 4))
Balanced_accuracy = (Sensitivity + Specificity) / 2

xlab = sprintf('Sensitivity = %0.3g     Specificity = %0.3g     Balanced Acc = %0.3g', Sensitivity, Specificity, Balanced_accuracy);
xlabel(xlab)
% Weighted_accuracy = 0.25*Sensitivity + 0.75*Specificity

% Accuracy = ( sum(toplot == 1) + sum(toplot == 3) ) / (sum(toplot == 1) + sum(toplot == 2) + sum(toplot == 3) + sum(toplot == 4) )



%%%% Display parameters
% Algo
disp(['Ground truth: ' compare_truth])
disp(['Automatic classifier: ' compare_with])

% Nb of IC
fprintf('Nb of IC per subject: %d (%s)\n',size(toplot,2),  nb_IC_of_interest)

% Which col of the csv
if isempty(which_col_in_csv)
    disp('Cardiac IC tested: ALL (sure + noisy)')
else
    disp(['Cardiac IC tested: ' which_col_in_csv])
end


%% Function to interact with the plot
function plotit(src,evt)

tmp = get(gcf,'UserData');
fs = tmp.fs;
compare_with = tmp.compare_with;
compare_truth = tmp.compare_truth;

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

switch compare_with
    case 'CARACAS'
        toprint = removefields(fs(ids).CARACAS.meas(icmp),{'Ampl_var'});
        fields = fieldnames(toprint);
        threshs_min = [0 0 0 0 0 0 2 5 0 0 35];
        threshs_max = [1/3 1/3 1/3 1/3 1/3 1/3 inf 100 1/3 1/3 90];
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
    case 'SASICARACAS'
        toprint = fs(ids).SASICARACAS.meas(icmp);
        cfg = fs(ids).SASICARACAS.cfg;
        fields = fieldnames(toprint);
        for i = 1:numel(fields)
            strtitle = [fields{i} ' = ' num2str(toprint.(fields{i}),2)];
            c = 'k';
            if isscalar(cfg.(['thresh_' fields{i}]))
                cfg.(['thresh_' fields{i}]) = [0 cfg.(['thresh_' fields{i}])];
            end
            c = ifelse(toprint.(fields{i}) < cfg.(['thresh_' fields{i}])(1),'r',c);
            c = ifelse(toprint.(fields{i}) > cfg.(['thresh_' fields{i}])(2),'r',c);
            x = xl(1)+(diff(xl)/3)*rem(i-1,round(numel(fields)/3)) +1;
            y = yl(1)-(floor((i)/(numel(fields)/3))+1)*diff(yl)/20-diff(yl)/10;
            text(x,y,strtitle,'HorizontalAlignment','left','VerticalAlignment','baseline','FontSize',8, 'Interpreter','none', 'color',c)
        end
    case 'ICLabel'
        classes = fs(ids).ICLabel.meas.classes;
        classifications = fs(ids).ICLabel.meas.classifications(icmp,:);
        for i = 1:numel(classes)
            strtitle = [classes{i} ' = ' num2str(classifications(i),2)];
            c = 'k';
            % Highlight the class with highest probability in red
            if classifications(i) == max(classifications)
                c = 'r';
            end
            x = xl(1)+(diff(xl)/3)*rem(i-1,round(numel(classes)/3)) +1;
            y = yl(1)-(floor((i)/(numel(classes)/3))+1)*diff(yl)/20-diff(yl)/10;
            text(x,y,strtitle,'HorizontalAlignment','left','VerticalAlignment','baseline','FontSize',8, 'Interpreter','none', 'color',c)
        end
        [~, max_idx] = max(classifications);
        title(['ICLabel: ' classes{max_idx} ' (' num2str(classifications(max_idx),2) ')'])
end

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

ids

corr = fs(ids).CORR.c(icmp)

% fs(ids).(compare_truth).meas(icmp)
% fs(ids).(compare_with).meas(icmp)

% cfg = [];
% cfg.viewmode  = 'component';
% cfg.layout    = layout;
% ft_databrowser(cfg, comp);

end
