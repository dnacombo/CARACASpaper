

% clear all
% close all

%% !% Create output directory for plots
setpath_ds003690

%% Path and load data
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

plot_dir = fullfile(dirout, '..', 'CardiClassif_SASICA', 'plots_compare');! USER CHOOSE PARAMETERS !!!!!

which_data = '_SASICA'; %''; %'_30comp';%   //

nb_IC_of_interest = 56;%'all'; %55;%45;% 'all' // 

which_col_in_csv = '_sure'; % '' // '_sure' // 'rej_noisy'

compare_truth = 'MANUAL';

% Define all parameter combinations to test
param_combinations = {'SASICARACAS_orig', 'SASICARACAS_absPT', 'SASICARACAS_abstemplate', 'SASICARACAS_NaNST', 'SASICARACAS_absPT_abstemplate', 'SASICARACAS_absPT_NaNST', 'SASICARACAS_abstemplate_NaNST', 'SASICARACAS_absPT_abstemplate_NaNST'};

% Also include standard comparisons
standard_comparisons = {'CARACAS', 'SASICARACAS', 'ICLabel', 'CORR'};
all_comparisons = [standard_comparisons, param_combinations];

% Create output directory for plots
setpath_ds003690
plot_dir = fullfile(dirout, '..', 'CardiClassif_SASICA', 'plots_compare');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% Loop through all comparisons
for i_comp = 1:length(all_comparisons)
    compare_with = all_comparisons{i_comp};
    
    fprintf('\n=== Processing comparison %d/%d: %s ===\n', i_comp, length(all_comparisons), compare_with);
    
    % Check if this comparison exists in the data
    if ~strcmp(compare_with, 'CARACAS') && ~strcmp(compare_with, 'ICLabel') && ~strcmp(compare_with, 'CORR')
        % Check if any file has this field
        has_field = false;
        for i_f = 1:min(5, numel(fs)) % Check first 5 files
            if isfield(fs(i_f), compare_with)
                has_field = true;
                break;
            end
        end
        if ~has_field
            fprintf('   Skipping %s - field not found in data\n', compare_with);
            continue;
        end
    end
    
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
    set(gcf,'UserData',struct('fs',fs, 'compare_with', compare_with, 'compare_truth',compare_truth, 'cfg_CARACAS', cfg_CARACAS))
    subplot(321)
    hi = imagesc(truthrej');
    title(compare_truth, 'interpreter', 'none')
    set(hi,"ButtonDownFcn",@plotit);
    subplot(322);
    hi = imagesc(withrej');
    title(compare_with, 'interpreter', 'none')
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
        toplot = toplot(:,1:nb_IC_of_interest);
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
    Sensitivity = sum(toplot == 1) / (sum(toplot== 1) + sum(toplot == 2));
    % Specificity = CR  / (CR + FA)
    Specificity = sum(toplot == 3) / (sum(toplot== 3) + sum(toplot == 4));
    Balanced_accuracy = (Sensitivity + Specificity) / 2;

    xlab = sprintf('Sensitivity = %0.3g     Specificity = %0.3g     Balanced Acc = %0.3g', Sensitivity, Specificity, Balanced_accuracy);
    xlabel(xlab)

    % Create overall title for the figure
    sgtitle(sprintf('Comparison: %s vs %s', compare_truth, compare_with), 'FontSize', 14, 'FontWeight', 'bold', 'interpreter', 'none');

    %%%% Display and save results
    fprintf('   Ground truth: %s\n', compare_truth);
    fprintf('   Automatic classifier: %s\n', compare_with);
    fprintf('   Sensitivity: %.3f\n', Sensitivity);
    fprintf('   Specificity: %.3f\n', Specificity);
    fprintf('   Balanced Accuracy: %.3f\n', Balanced_accuracy);

    % Maximize figure before saving
    set(gcf, 'WindowState', 'maximized');
    
    % Save figure
    fig_name = sprintf('SDT_comparison_%s_vs_%s', compare_truth, compare_with);
    fig_path = fullfile(plot_dir, [fig_name '.fig']);
    saveas(gcf, fig_path);
    fig_path = strrep(fig_path,'.fig','.png');
    saveas(gcf, fig_path);
    fprintf('   Figure saved: %s\n', fig_path);
    
    % Save performance metrics
    results(i_comp).compare_with = compare_with;
    results(i_comp).compare_truth = compare_truth;
    results(i_comp).sensitivity = Sensitivity;
    results(i_comp).specificity = Specificity;
    results(i_comp).balanced_accuracy = Balanced_accuracy;
    results(i_comp).n_subjects = size(toplot, 1);
    results(i_comp).n_components = size(toplot, 2);
end

%% Save summary results
summary_file = fullfile(plot_dir, 'comparison_summary.mat');
save(summary_file, 'results');
fprintf('\nSummary results saved: %s\n', summary_file);

% Display summary table
fprintf('\n=== SUMMARY OF ALL COMPARISONS ===\n');
fprintf('%-35s | Sens. | Spec. | Bal.Acc\n', 'Method');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:length(results)
    fprintf('%-35s | %.3f | %.3f | %.3f\n', ...
        results(i).compare_with, ...
        results(i).sensitivity, ...
        results(i).specificity, ...
        results(i).balanced_accuracy);
end


%% Function to interact with the plot
function plotit(src,evt)

tmp = get(gcf,'UserData');
fs = tmp.fs;
compare_with = tmp.compare_with;
compare_truth = tmp.compare_truth;
cfg_CARACAS = tmp.cfg_CARACAS;

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

title(['comp' num2str(icmp) ' ' fs(ids).sub,' ', fs(ids).task,' ', fs(ids).run], 'fontsize',16, 'Interpreter', 'none');


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
        
        % Use CARACAS_rethresh function to determine which criteria failed
        % Create a temporary structure with only the current component
        temp_CARACAS = fs(ids).CARACAS;
        temp_CARACAS.meas = fs(ids).CARACAS.meas(icmp);
        dummy_rej = zeros(1, 1);
        [~, failed_criteria] = CARACAS_rethresh(dummy_rej, temp_CARACAS, cfg_CARACAS);
        
        for i = 1:numel(fields)
            strtitle = [fields{i} ' = ' num2str(toprint.(fields{i}),2)];
            c = 'k';
            
            % Color red if this specific field failed the criteria
            if isfield(failed_criteria, fields{i}) && failed_criteria.(fields{i})
                c = 'r';
            end
            
            % Calculate position for clean column layout
            ncols = min(3, numel(fields)); % Maximum 3 columns
            nrows = ceil(numel(fields) / ncols);
            col = rem(i-1, ncols) + 1;
            row = ceil(i / ncols);
            
            x = xl(1) + (col-1) * diff(xl) / ncols + diff(xl) / (ncols * 20);
            y = yl(1) - diff(yl) * 0.1 - (row-1) * diff(yl) / 10;
            
            text(x, y, strtitle, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8, 'Interpreter', 'none', 'color', c)
        end
        title([ifelse(fs(ids).CARACAS.rej(icmp),'CARACAS ',''),ifelse(fs(ids).MANUAL.rej(icmp), ' MANUAL ' ,'')], 'interpreter', 'none')
    case {'SASICARACAS', 'nuSASICARACAS'}
        toprint = fs(ids).SASICARACAS.meas(icmp);
        cfg = fs(ids).SASICARACAS.cfg;
        fields = fieldnames(toprint);
        
        % Use CARACAS_rethresh function to determine which criteria failed
        % Create a temporary structure with only the current component
        temp_SASICARACAS = fs(ids).SASICARACAS;
        temp_SASICARACAS.meas = fs(ids).SASICARACAS.meas(icmp);
        dummy_rej = zeros(1, 1);
        [~, failed_criteria] = CARACAS_rethresh(dummy_rej, temp_SASICARACAS, cfg);
        
        for i = 1:numel(fields)
            strtitle = [fields{i} ' = ' num2str(toprint.(fields{i}),2)];
            c = 'k';
            
            % Color red if this specific field failed the criteria
            if isfield(failed_criteria, fields{i}) && failed_criteria.(fields{i})
                c = 'r';
            end
            
            % Calculate position for clean column layout
            ncols = min(3, numel(fields)); % Maximum 3 columns
            nrows = ceil(numel(fields) / ncols);
            col = rem(i-1, ncols) + 1;
            row = ceil(i / ncols);
            
            x = xl(1) + (col-1) * diff(xl) / ncols + diff(xl) / (ncols * 20);
            y = yl(1) - diff(yl) * 0.1 - (row-1) * diff(yl) / 10;
            
            text(x, y, strtitle, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8, 'Interpreter', 'none', 'color', c)
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
            
            % Calculate position for clean column layout
            ncols = min(3, numel(classes)); % Maximum 3 columns
            nrows = ceil(numel(classes) / ncols);
            col = rem(i-1, ncols) + 1;
            row = ceil(i / ncols);
            
            x = xl(1) + (col-1) * diff(xl) / ncols + diff(xl) / (ncols * 20);
            y = yl(1) - diff(yl) * 0.05 - (row-1) * diff(yl) / (nrows + 5);
            
            text(x, y, strtitle, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 8, 'Interpreter', 'none', 'color', c)
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
