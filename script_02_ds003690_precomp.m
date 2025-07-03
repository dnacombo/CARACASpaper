function script_02_ds003690_precomp(i_f)

if not(exist('i_f','var'))
    i_f = 1:375;
end
plot_epochs = 0;
plot_psd = 0;
plot_comp = 0;
force_recomp = 0;

corthresh = .6; % correlation threshold component time course with ECG

which_data = '_SASICA'; % '' // '_30comp'

setpath_ds003690
ft_warning('off','FieldTrip:dataContainsNaN')


if not(exist(fullfile(dirout,'AllFilesAndScoresList.mat'), 'file'))
    fs = flister('sub.*/(?<sub>sub-[^_]+)_(?<task>task-[^_]+)_(?<run>run-[^_]+)_(?<mod>eeg).set','dir',dirbids);
    save(fullfile(dirout,'AllFilesAndScoresList.mat'),'fs')
else
    load(fullfile(dirout,'AllFilesAndScoresList.mat'), 'fs')
end

% fs = flister('meg.fif', 'dir',fullfile(dirroot), 'recurse',0);
% fs = fs(1);

% remove sub 1 with strange cardiac rhytms
% fs = flist_select(fs,'sub', 'sub-AB10', 'inv');
%% 
for i_f = i_f%1:numel(fs)
    fprintf('#############################################\n')
    fprintf('######### %s_%s_%s #########\n',fs(i_f).sub,fs(i_f).task,fs(i_f).run)
    fprintf('#############################################\n')
    %% read and segment data
    cfg = [];
    cfg.dataset = fs(i_f).name;
    this_outdir = fullfile(dirout,fs(i_f).sub, fs(i_f).mod);
    mymkdir(this_outdir);
    lock = fullfile(this_outdir,'lock');
    if exist(lock,'file')
        continue
    else
        % fclose(fopen(lock,'a'));
    end

    this_file_eeg = fullfile(this_outdir,[myfileparts(fs(i_f).name,'f') '.mat']);
    this_file_comp = strrep(this_file_eeg,'_eeg','_comp');
    fs(i_f).compf = this_file_comp;
    fs(i_f).eegf = this_file_eeg;

    try
        if force_recomp
            error('planned force recomp')
        end
        data = load(this_file_eeg,'-mat');
        comp = load(this_file_comp,'-mat');
    catch
        data = [];
        comp = [];
    end

    if isempty(comp)
        % average reference
        cfg.reref = 'yes';
        cfg.channel =  {'all' '-M1' '-M2'};
        cfg.refmethod = 'avg';
        cfg.refchannel = 'all';%{'M1', 'M2'};
        % Filter
        cfg.hpfilter = 'yes';
        cfg.hpfreq = .1;
        cfg.hpfilttype = 'firws';
        cfg.lpfilter = 'yes';
        cfg.lpfilter = 35;
        cfg.lpfilttype = 'firws';

        data = ft_preprocessing(cfg);

        cfg.trialdef.eventtype = 'cue';
        cfg.trialdef.prestim = 0.2;
        cfg.trialdef.poststim = 6;
        cfg = ft_definetrial(cfg);

        cfg.demean = 'yes';
        cfg.baselinewindow = [-.2 0];
        data = ft_redefinetrial(cfg,data);
        %%
        if plot_epochs
            cfg = [];
            ft_databrowser(cfg,data);
        end

        %%
        cfg = [];
        cfg.channel = 'eeg';
        data.elec.coordsys = 'EEGLAB';

        layout = ft_prepare_layout(cfg,data);

        %%
        cfg = [];
        cfg.channel = 'EEG';
        % cfg.numcomponent = 30;

        rng(123);
        comp = ft_componentanalysis(cfg,data);
        save(this_file_comp,'-fromstruct',comp);
        save(this_file_eeg,'-fromstruct',data);
    end


    %% control plot psd
    if plot_psd
        cfg = [];
        cfg.output = 'pow';
        cfg.method = 'mtmfft';
        cfg.taper   = 'boxcar';
        cfg.channel = 'eeg';
        cfg.foi     = 0.5:1:45;
        base_freq   = ft_freqanalysis(cfg, data);

        figure(2);clf;
        semilogy(base_freq.freq,base_freq.powspctrm)
        hold on
        m = mean(base_freq.powspctrm,1);
        sd = std(base_freq.powspctrm,[],1);
        err = [m+sd;m-sd];
        yl = ylim;
        err(err<0) = yl(1);
        semilogy(base_freq.freq,m, 'k', 'LineWidth',2)
        semilogy(base_freq.freq,err, 'k:')
        xlabel('Freq')
        ylabel('Pow')
    end
    %% plot components
    if plot_comp
        %%
        cfg = [];
        cfg.layout = layout;
        cfg.component = 1:20;
        cfg.channel = 'eeg';
        ft_topoplotIC(cfg,comp)

        %%
        cfg = [];
        cfg.viewmode  = 'component';
        cfg.layout    = layout;
        ft_databrowser(cfg, comp);

    end
    cfgtmp = [];
    if isfield(comp, 'sampleinfo')
        cfgtmp.trl = [1, comp.sampleinfo(end,2), 0];
    else
        cfgtmp.trl = [1, sum(cellfun(@numel,comp.time)), 0]; % Start, End, Offset
    end
    comp_continu = ft_redefinetrial(cfgtmp, comp);

    %% CARACAS in SASICA
    cfg = [];
    cfg.channel = 'eeg';
    EEG = comp2eeglab(cfg, comp, ft_selectdata(cfg,data));

    cfg_SASICA = [];
    cfg_SASICA.CARACAS.enable = 1;
    cfg_SASICA.opts.noplot = 1;
    cfg_SASICA.opts.noplotselectcomps = 1;

    SASICCARACAS = eeg_SASICA(EEG,cfg_SASICA);
    fs(i_f).SASICARACAS.cfg = SASICCARACAS.reject.SASICA.icaCARACAS_cfg;
    fs(i_f).SASICARACAS.meas = SASICCARACAS.reject.SASICA.icaCARACAS;
    fs(i_f).SASICARACAS.rej = SASICCARACAS.reject.gcompreject;
    
    %% CARACAS

    cfg = [];
    cfg.method_chosen = 'absolute_threshold';   %'absolute_threshold' (method 1) or 'mean_std' (method 2)
    cfg.plot_heart_IC = 0;                      % 1 or 0 (to plot the IC labelled as cardiac)
    cfg.path_output = this_outdir; % path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
    cfg.file_info = myfileparts(fs(i_f).name,'f');% name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)
    cfg.nb_IC_wanted = 3;                       % number of IC selected for each metric (kurtosis, skewness...) [default: 3, to select the top 3 IC for each metric]
    cfg.bpm_min = 45;                           % expected heart beat per min, for sanity check [default: 45 and 90]
    cfg.bpm_max = 90;
    cfg.threshold_cond_IC_method1 = .5;         % minimum proportion of conditions that must be met in order that an IC could be considered as a potential heart IC [default: 0.5, so if method_chosen == 'absolute_threshold', an IC must be in the top 3 for at least 50% of the metrics]
    cfg.threshold_std_method2 = 2.5;            % if method_chosen == 'mean_std', an IC will be considered as a potential heart IC if its proportion of conditions met (i.e., its score) is above mean(all_score) + threshold_std_method2 * std(all_score) [default: 2.5]
    cfg.min_recording_duration_sec = 20;        % minimum duration (in sec) of the IC timecourse (default: 20]
    cfg.mini_bouts_duration_for_SignalAmplRange = 10; % for sanity check (avoids false positive): the time course of a potential heart IC must be ~regular. The timecourse will be divided into mini-segments of this duration, and we will check that the amplitude between these mini-bouts is ~similar. [default: 10]
    cfg.threshold_regularity_signal_minmax = 1.5; % For each mini-bout, the averaged signal amplitude is computed. The IC timecourse will be considered as irregular if: (max(Mean_Amp_minibout) - min(Mean_Amp_minibout)) / min(Mean_Amp_minibout) > threshold_regularity_signal_minmax [default: 1.5]

    [tmp, meas] = CARACAS(cfg, comp_continu);
    fs(i_f).CARACAS.meas = meas;
    fs(i_f).CARACAS.rej = false(1,numel(comp.label));
    fs(i_f).CARACAS.rej(tmp) = 1;


    %% IClabel
    activate_matconvnet();
    EEG = pop_iclabel(EEG,'default');
    ICLabel_results = EEG.etc.ic_classification.ICLabel;
    fs(i_f).ICLabel.meas = struct();
    fs(i_f).ICLabel.meas.classes = ICLabel_results.classes;
    fs(i_f).ICLabel.meas.classifications = ICLabel_results.classifications;
    cls = fs(i_f).ICLabel.meas.classifications(:,~strcmp(fs(i_f).ICLabel.meas.classes,'Heart'));
    m = max(cls,[],2);
    tmp = fs(i_f).ICLabel.meas.classifications(:,strcmp(fs(i_f).ICLabel.meas.classes,'Heart'));
    fs(i_f).ICLabel.rej = tmp' > m';


    %% correlation with EKG

    EKGchan = chnb('ekg',data.label);

    EKG = cat(2,data.trial{:});
    EKG = EKG(EKGchan,:,:)';
    ICs = cat(2,comp.trial{:})';

    c  = abs(corr(ICs,EKG))';
    rej = c > corthresh ;

    fs(i_f).CORR.c = c;
    fs(i_f).CORR.rej = rej;

    %% SDT
    fs(i_f).pC = sum(fs(i_f).CORR.rej == fs(i_f).CARACAS.rej) / numel(fs(i_f).CORR.rej);
    fs(i_f).dp = PAL_SDT_MAFC_PCtoDP(fs(i_f).pC,numel(fs(i_f).CORR.rej));
    fs(i_f).H = sum(fs(i_f).CORR.rej & fs(i_f).CARACAS.rej);
    fs(i_f).FA = sum(~ fs(i_f).CORR.rej & fs(i_f).CARACAS.rej);
    
    % delete(lock)
end
f = fs(i_f);
save(fullfile(dirout,f.sub,sprintf('%s_%s_%s_FilesAndScoresList.mat',f.sub, f.task, f.run)),'f')

% %% ROC analysis
% [X, Y, T, AUC] = perfcurve(double([fs.CORR]), double([fs.CARACAS]), 1);
% figure(587934);clf
% plot(X, Y, 'Marker','+');
% xlabel('False Positive')
% ylabel('True Positive')
% title(['Classification perf. AUC = ' num2str(AUC)])
