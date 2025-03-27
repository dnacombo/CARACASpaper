
dirroot = '/network/iss/cenir/analyse/meeg/CARACAS/Test_Max/';
dircode = fullfile(dirroot,'code');
dirbids       = fullfile(dirroot,'ds002718');
dirderiv = fullfile(dirbids,'derivatives');
dirout = fullfile(dirderiv,'CardiClassif');

addpath(fullfile(dircode,'MiscMatlab'))
addpath(fullfile(dircode,'MiscMatlab', 'stats'))
addpath(fullfile(dircode,'SASICA'))
rm_frompath('eeglab') % avoid interference with other versions
rm_frompath('fieldtrip')
% addpath(fullfile(dircode,'MiscMatlab', 'eeglab'))
% addpath(fullfile(dircode,'eeglab'))
addpath(fullfile(dircode,'fieldtrip'))
addpath(fullfile(dircode,'fieldtrip', 'external','eeglab'))
ft_defaults

% start EEGLAB
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

addpath(fullfile(dircode,'Cardiac_IC_labelling'));

cfg_SASICA = SASICA('getdefs');
addpath(genpath(fullfile(dircode, "SASICA",'eeglab')))


fs = flister('derivatives/eeglab/sub.*/(?<sub>sub-\d+)_(?<task>task-\w+)_(?<step>02-ICA)_(?<mod>eeg).set','dir',dirderiv);
% fs = fs(1);
i_f = 1;
%%
heart_IC = [];
%% CARACAS
for i_f = 1:numel(fs)
    %%
    EEG = load(fs(i_f).name, '-mat');
    %%
    comp = eeglab2fieldtrip(EEG,'comp');
    comp.fsample = EEG.srate;

    cfg = [];
    cfg.method_chosen = 'absolute_threshold';   %'absolute_threshold' (method 1) or 'mean_std' (method 2)
    cfg.plot_heart_IC = 0;                      % 1 or 0 (to plot the IC labelled as cardiac)
    cfg.path_output = fullfile(dirout,fs(i_f).sub, fs(i_f).mod); % path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
    cfg.file_info = EEG.setname;                % name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)
    cfg.nb_IC_wanted = 3;                       % number of IC selected for each metric (kurtosis, skewness...) [default: 3, to select the top 3 IC for each metric]
    cfg.bpm_min = 45;                           % expected heart beat per min, for sanity check [default: 45 and 90]
    cfg.bpm_max = 90;
    cfg.threshold_cond_IC_method1 = .5;         % minimum proportion of conditions that must be met in order that an IC could be considered as a potential heart IC [default: 0.5, so if method_chosen == 'absolute_threshold', an IC must be in the top 3 for at least 50% of the metrics]
    cfg.threshold_std_method2 = 2.5;            % if method_chosen == 'mean_std', an IC will be considered as a potential heart IC if its proportion of conditions met (i.e., its score) is above mean(all_score) + threshold_std_method2 * std(all_score) [default: 2.5]
    cfg.min_recording_duration_sec = 20;        % minimum duration (in sec) of the IC timecourse (default: 20]
    cfg.mini_bouts_duration_for_SignalAmplRange = 10; % for sanity check (avoids false positive): the time course of a potential heart IC must be ~regular. The timecourse will be divided into mini-segments of this duration, and we will check that the amplitude between these mini-bouts is ~similar. [default: 10]
    cfg.threshold_regularity_signal_minmax = 1.5; % For each mini-bout, the averaged signal amplitude is computed. The IC timecourse will be considered as irregular if: (max(Mean_Amp_minibout) - min(Mean_Amp_minibout)) / min(Mean_Amp_minibout) > threshold_regularity_signal_minmax [default: 1.5]

    mymkdir(cfg.path_output);
    [tmp] = CARACAS(cfg, comp);
    heart_IC(i_f).CARACAS = false(size(comp.label));
    heart_IC(i_f).CARACAS(tmp) = 1;
end
%% correlation with EKG
for i_f = 1:numel(fs)
    %%
    EEG = load(fs(i_f).name, '-mat');
    
    EKGchan = chnb('ekg',{EEG.chanlocs.type});
    cfg_SASICA.chancorr.enable = 1;
    cfg_SASICA.chancorr.channames = EKGchan;
    cfg_SASICA.chancorr.corthresh = .8;
    cfg_SASICA.opts.noplot = 1;
    cfg_SASICA.opts.noplotselectcomps = 1;

    EEG = eeg_SASICA(EEG, cfg_SASICA);
    heart_IC(i_f).SASICA = EEG.reject.gcompreject;

end
% %% CARACAS in SASICA
% for i_f = 1:numel(fs)
%     %%
%     EEG = load(fs(i_f).name, '-mat');
% 
%     cfg_SASICA.CARACAS.enable = 1;
%     cfg_SASICA.opts.noplot = 1;
%     cfg_SASICA.opts.noplotselectcomps = 1;
% 
%     EEG = eeg_SASICA(EEG, cfg_SASICA);
%     heart_IC(i_f).CARACASICA = EEG.reject.gcompreject;
% 
% end
%%
[X, Y, T, AUC] = perfcurve(double([heart_IC.SASICA]), double([heart_IC.CARACAS]), 1);
figure(587934);clf
plot(X, Y);
xlabel('False Positive')
ylabel('True Positive')
title(['Classification perf. AUC = ' num2str(AUC)])