function script_02_ds003690_precomp(i_f)

if not(exist('i_f','var'))
    i_f = 1;
end
plot_epochs     = 0;
plot_psd        = 0;
plot_comp       = 0;
force_recomp    = 0;

quick_update        = 0;
update_CARACAS      = 0;
update_SASICARACAS  = 0;
update_nuSASICARACAS= 1;
update_ICLabel      = 0;
update_CORR         = 0;

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
% fs = flist_select(fs,'sub', 'sub-AB46', 'task', 'task-simpleRT', 'run', 'run-1');
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

    if ~quick_update
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

    end

    if update_SASICARACAS
        %% CARACAS in SASICA
        cfg = [];
        cfg.channel = 'eeg';
        EEG = comp2eeglab(cfg, comp, ft_selectdata(cfg,data));

        cfg_SASICA = SASICA('getdefs');
        cfg_SASICA.CARACAS.enable = 1;
        cfg_SASICA.opts.noplot = 1;
        cfg_SASICA.opts.noplotselectcomps = 1;

        SASICCARACAS = eeg_SASICA(EEG,cfg_SASICA);
        fs(i_f).SASICARACAS.cfg = SASICCARACAS.reject.SASICA.icaCARACAS_cfg;
        fs(i_f).SASICARACAS.meas = SASICCARACAS.reject.SASICA.icaCARACAS;
        fs(i_f).SASICARACAS.rej = SASICCARACAS.reject.gcompreject;
    end
    if update_nuSASICARACAS
        %% CARACAS in SASICA after modifying eeg_SASICA to use option absPT for heart_peak_detect
        cfg = [];
        cfg.channel = 'eeg';
        EEG = comp2eeglab(cfg, comp, ft_selectdata(cfg,data));

        cfg_SASICA = SASICA('getdefs');
        cfg_SASICA.CARACAS.enable = 1;
        cfg_SASICA.opts.noplot = 1;
        cfg_SASICA.opts.noplotselectcomps = 1;

        % Test all 8 combinations of the three parameters
        param_combinations = [
            0, 0, 0;  % absPT=0, abstemplate=0, NaNST=0 (original)
            1, 0, 0;  % absPT=1, abstemplate=0, NaNST=0
            0, 1, 0;  % absPT=0, abstemplate=1, NaNST=0
            0, 0, 1;  % absPT=0, abstemplate=0, NaNST=1
            1, 1, 0;  % absPT=1, abstemplate=1, NaNST=0
            1, 0, 1;  % absPT=1, abstemplate=0, NaNST=1
            0, 1, 1;  % absPT=0, abstemplate=1, NaNST=1
            1, 1, 1;  % absPT=1, abstemplate=1, NaNST=1
        ];
        
        param_names = {'orig', 'absPT', 'abstemplate', 'NaNST', 'absPT_abstemplate', 'absPT_NaNST', 'abstemplate_NaNST', 'absPT_abstemplate_NaNST'};
        
        for i_comb = 1:size(param_combinations, 1)
            % Set parameters for this combination
            cfg_SASICA.CARACAS.cfg_peak.absPT = param_combinations(i_comb, 1);
            cfg_SASICA.CARACAS.cfg_peak.abstemplate = param_combinations(i_comb, 2);
            cfg_SASICA.CARACAS.cfg_peak.NaNST = param_combinations(i_comb, 3);
            
            % Run SASICA with current parameter combination
            SASICCARACAS = eeg_SASICA(EEG, cfg_SASICA);
            
            % Store results with descriptive field name
            field_name = ['SASICARACAS_' param_names{i_comb}];
            fs(i_f).(field_name).cfg = SASICCARACAS.reject.SASICA.icaCARACAS_cfg;
            fs(i_f).(field_name).meas = SASICCARACAS.reject.SASICA.icaCARACAS;
            fs(i_f).(field_name).rej = SASICCARACAS.reject.gcompreject;
        end
    end
    if update_CARACAS
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
    end

    if update_ICLabel
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
    end

    if update_CORR
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
    end

    % delete(lock)
    f = fs(i_f);
    save(fullfile(dirout,f.sub,sprintf('%s_%s_%s_FilesAndScoresList.mat',f.sub, f.task, f.run)),'f')

end
% %% ROC analysis
% [X, Y, T, AUC] = perfcurve(double([fs.CORR]), double([fs.CARACAS]), 1);
% figure(587934);clf
% plot(X, Y, 'Marker','+');
% xlabel('False Positive')
% ylabel('True Positive')
% title(['Classification perf. AUC = ' num2str(AUC)])
end

function [EEG,cfg] = comp2eeglab(cfg,comp,data)

% create an EEG structure based on comp.

EEG = eeg_emptyset;
EEG.setname = 'internal';
EEG.nbchan = numel(comp.topolabel);
if EEG.nbchan == 0
    error('No more channels here.')
end
EEG.trials = numel(comp.trial);
EEG.pnts = size(comp.trial{1},2);
EEG.srate = comp.fsample;
EEG.xmin = comp.time{1}(1);
EEG.xmax = comp.time{end}(end);
EEG.times = comp.time{1}*1000;
if isscalar(unique(cellfun(@(x) size(x,2),comp.trial)))
    EEG.icaact = cat(3,comp.trial{:});
else
    warning('Trials have unequal length. Catenating for display')
    EEG.icaact = cat(2,comp.trial{:});
end
EEG.icawinv = comp.topo;
EEG.icaweights = pinv(EEG.icawinv);
EEG.icasphere  = eye(size(EEG.icaweights,2));


EEG.chanlocs = struct();
if isfield(cfg,'layout')
    cfg.layout = ft_prepare_layout(cfg);
else
    if not(isfield(data,'elec'))
        warning('No layout provided. Topographies may be inaccurate')
    end
end
for i = 1:EEG.nbchan
    EEG.chanlocs(i).labels = comp.topolabel{i};
    % attempt to create a chanlocs
    if isfield(cfg,'layout')
        ichan = chnb(comp.topolabel{i},cfg.layout.label);
        if ~isempty(ichan)
            [EEG.chanlocs(i).X] = cfg.layout.pos(ichan,1);
            [EEG.chanlocs(i).Y] = cfg.layout.pos(ichan,2);
            [EEG.chanlocs(i).Z] = 1;
        end
    elseif exist('data','var') && isfield(data,'elec')
        ichan = chnb(comp.topolabel{i},data.elec.label);
        if isfield(data.elec,'pnt')
            EEG.chanlocs(i).X = data.elec.pnt(ichan,1);
            EEG.chanlocs(i).Y = data.elec.pnt(ichan,2);
            EEG.chanlocs(i).Z = data.elec.pnt(ichan,3);
        elseif isfield(data.elec,'chanpos')
            EEG.chanlocs(i).X = data.elec.chanpos(ichan,1);
            EEG.chanlocs(i).Y = data.elec.chanpos(ichan,2);
            EEG.chanlocs(i).Z = data.elec.chanpos(ichan,3);
        end            
    end
end
EEG.chanlocs = convertlocs(EEG.chanlocs,'cart2all');
mr = max([EEG.chanlocs.radius]);
for i = 1:numel(EEG.chanlocs)
    EEG.chanlocs(i).radius = EEG.chanlocs(i).radius * .5/ mr;
end
EEG.chaninfo.nosedir = '+Y';

if not(exist('data','var'))
    EEG.data = reshape(EEG.icawinv * EEG.icaact(:,:),EEG.nbchan,EEG.pnts,EEG.trials);
else
    if isscalar(unique(cellfun(@(x) size(x,2),comp.trial)))
        EEG.data = cat(3,data.trial{:});
    else
        EEG.data = cat(2,data.trial{:});
        EEG.trials = 1;
    end
end
EEG.pnts = size(EEG.data,2);
EEG.icachansind = 1:size(EEG.data,1);

if isfield(cfg,'reject')
    EEG.reject = cfg.reject;
end
end