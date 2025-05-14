%% DOCSTRING

% Pierre Champetier (November 2024)


%% Ressources

% EEGLAB tutorial for ICA:
% --> https://eeglab.org/tutorials/06_RejectArtifacts/RunICA.html#inspecting-ica-components

% SASICA (creator = Maximilien Chaumon, ICM):
% -->https://github.com/dnacombo/SASICA


%% Clear Workspace

clear all;
close all;


%% CONSTANT PARAMETERS

% Save .png of ICA quality check (IC toppographies + IC timecourse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_png = 'n';

% Which data
%%%%%%%%%%%%
% if ~isempty(findstr(pwd,'lustre')) || ~isempty(findstr(pwd,'iss02'))
    which_data = 'RS';
% else
%     prompt = {"Which data you want to pre-process? (RS for resting-state, T for FCSRT task, RST for both, RSTW for RS+WordList+FCSRT)" };
%     dlgtitle = 'Data';
%     dims = [1 110];
%     definput = {'RS'};
%     answer = inputdlg(prompt,dlgtitle,dims,definput);
%     which_data = char(answer)
% end



% Scalp EEG channels
%%%%%%%%%%%%%%%%%%%%
% Do you want to keep prefrontal EEG channels? (k for keeping all, r for removing part of them, w for without all prefrontal channels, a for keeping all channels including neck and cheeks)
prefrontal_EEG = 'keep_prefrontal_EEG'; % 'all_EEG_channels' / 'keep_prefrontal_EEG' / 'remove_part_prefrontal_EEG', 'without_prefrontal_EEG'


%% Paths
dirroot = '/network/iss/cenir/analyse/meeg/CARACAS/Test_Max/';
dircode = fullfile(dirroot,'code');
dirbids       = fullfile(dirroot,'../Test_Pierre/Data');
dirderiv = fullfile(dirbids,'derivatives');
dirout = fullfile(dirderiv,'CardiClassif');

addpath(fullfile(dircode,'MiscMatlab'))
addpath(fullfile(dircode,'MiscMatlab', 'stats'))
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


%% Define scalp electrodes

% run z_electrodes_toremove_EGI_PC.m;

%% FOR LOOP BEGINNING

% Define sessions
list_session = {'M0', 'M12', 'M24', 'M36', 'M48', 'M60'};

% Define runs
if strcmp(which_data, 'RS')
    list_run = [1,7];
elseif strcmp(which_data, 'T')
    list_run = [3:6];
elseif strcmp(which_data, 'RST')
    list_run = [1,3:7];
elseif strcmp(which_data, 'W')
    list_run = [2];
elseif strcmp(which_data, 'RSTW')
    list_run = [1:7];
end

% Initialize variables
sujets_erreurs = [];


%%%%%%%%%%%%%%%%%%
table_for_CHECKING = [];
%%%%%%%%%%%%%%%%%%

matrice_all = [];
zscore_ALL_heart_IC = cell(1,2);
zscore_ALL_heart_IC_sort = cell(1,2);
FINAL_ALL_FEATURES_table= [];
recording = 0;
classification_heart_for_ROC = [];
IC_ori = [];
sujet = 1;
%%
for sujet = 1:370
%%
    % ########################    $$$$$$$$$$$$******************************$$$$$$$$$$$$$$$$$$$$$$$$$$$$     ###################################################
    % if sujet == 14 || sujet == 13 || sujet == 20
    %     continue
    % end
    % ########################    $$$$$$$$$$$$******************************$$$$$$$$$$$$$$$$$$$$$$$$$$$$     ###################################################

    % Create complete name of the subject and display it
        sujet_name = sprintf('%03d',sujet);
    


    for session = 1:length(list_session)
        session_iter = char(list_session(session));

        for run_ = 1:length(list_run)
            run_iter = list_run(run_);

            %%%%%%%%%%%%%%%%%% CHOOSE A SPECIFIC SUBJECT FOR DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % list_sujet = [5];
            % if ismember(sujet, list_sujet) == 0
            %     continue
            % end
            % if strcmp(session_iter, 'M0') == 0
            %     continue
            % end
            % if run_ ~= 1
            %     continue
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if run_iter == 1 || run_iter == 7
                task_iter = 'RS';
            elseif run_iter == 2
                task_iter = 'WordList';
            else
                task_iter = 'FCSRT';
            end

            file_info = ['S' sujet_name '_' session_iter '_run_' num2str(run_iter)];
            display(['/// Analyzing ' file_info ' \\\'])

            % Check if subject has been preprocessed
            filename_output = [dirout '/ICA_sub-' sujet_name '_ses-' session_iter '_' task_iter '-run' num2str(run_iter) '.mat'];

            % List of EEG files (inputs)
            files_EEG = dir(dirbids);

            %% I - Load EEG data

            % 1) Find EEG file
            good_EEG = [];
            for i = 1:size(files_EEG,1)
                if contains(files_EEG(i,1).name, strcat('ICA_sub-', sujet_name, '_ses-', session_iter, '_', task_iter, '-run', num2str(run_iter), '')) == 1
                    good_EEG = [good_EEG, i];
                end
            end


            % 2) Check that only 1 EEG file has been found
            if size(good_EEG,2) == 0
                sujets_erreurs = [sujets_erreurs sujet];
                continue
            elseif size(good_EEG,2) > 1
                sujets_erreurs = [sujets_erreurs sujet];
                continue
            end


            % 3) If only 1 EEG file, load it
            load(strcat(dirbids, filesep, files_EEG(good_EEG,1).name));



            % 4) Update recording number
            recording = recording + 1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% V - ICA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      


            %% V - 3) Visual inspection of the topography and time course of the IC

            %%%%% To display epoch 2 (if RS) because subjects have eyes closed in epoch 1 during RS so no blinks (to compare with IC time course)
            if strcmp(task_iter, 'RS')
                epoch_to_display = 2;
            else
                epoch_to_display = 1;
            end


            if strcmp(save_png, 'y')

                % 1) Look at IC topography + timecourse
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cfg          = [];
                cfg.channel  = [1:10]; % components to be plotted
                cfg.viewmode = 'component';
                cfg.layout   = layout_short;

                cfg.trl(:, 1) = data_before_ICA.sampleinfo(epoch_to_display, 1);  % Start of each trial
                cfg.trl(:, 2) = data_before_ICA.sampleinfo(epoch_to_display, 2);  % End of each trial
                cfg.trl(:, 3) = 0;  % Offset

                % If only 1 epoch (for FCSRT and WordList), plot only the first 30s
                if (cfg.trl(:, 2) - cfg.trl(:, 1) ) / data_before_ICA.fsample > 30
                    cfg.trl(:, 2) = cfg.trl(:, 1) + 30*data_before_ICA.fsample;
                end

                % Plot
                try
                    ft_databrowser(cfg, comp)
                catch
                    continue
                end
                %%%%%%%%%%% SAVE fig as .png %%%%%%%%%%%
                filename = strcat(dirout, filesep, file_info, '_A_TimeCourse.png');
                % saveas(gcf, filename);
                % close(gcf);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                % 2) Look only at IC topography
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cfg          = [];
                cfg.component = 1:10;
                cfg.comment   = 'no';
                cfg.layout   = layout_short;
                ft_topoplotIC(cfg, comp)

                %%%%%%%%%%% SAVE fig as .png %%%%%%%%%%%
                filename = strcat(dirout, filesep, file_info, '_A_Topo.png'); % ADD BY PIERRE
                saveas(gcf, filename);
                close(gcf);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end


            %% V - 4) Identify heart IC (using EEGLAB)

            method_reject_cardiac_IC_parameters = [];
            method_reject_cardiac_IC_parameters.pop_runica_ica_label = [];
            method_reject_cardiac_IC_parameters.Pierre_fct = [];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % A- METHOD PIERRE (based on IC timecourse only)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            cfg = [];
            cfg.method_chosen = 'absolute_threshold'; %'absolute_threshold' (method 1) or 'mean_std' (method 2)
            cfg.plot_heart_IC = 0;
            cfg.path_output = dirout;
            cfg.file_info = file_info;
            

            % NEW FORMAT 2025
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [aaa_parameters_find_heart_IC, output_for_zscore_corMatrix_ROC, output_for_user] = CARACAS(cfg, comp);

            % Extract outputs
            heart_IC = output_for_user.heart_IC;

            z_score_heart_IC = output_for_zscore_corMatrix_ROC.z_score_heart_IC;
            correlationMatrix = output_for_zscore_corMatrix_ROC.correlationMatrix;
            FINAL_ALL_FEATURES_table = output_for_zscore_corMatrix_ROC.FINAL_ALL_FEATURES_table;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TO KNOW THE COND VERIFIED AMONG THE 12 COND
                % idx_cond_ok = find(table2array(table_IC_cond(rejected_comps_heart_PIERRE,:)) == 1);
                % metric_cond_ok = table_IC_cond.Properties.VariableNames(idx_cond_ok);
                % % Remove IC, bpm and SignalAmplRange name
                % if isequal(metric_cond_ok(1), {'IC'})
                %     metric_cond_ok = metric_cond_ok(2:end);
                % end
                % metric_cond_ok = metric_cond_ok(1:end-2);
        % % % % % % % % % % % % % % % % % % % % % 



        
        % Extract the classification result (0=non-cardiac, 1=cardiac) for all IC
        % % % % % % % % % % % % % % % % 
        classification_iter = repmat(0,1,length(comp.label));
        if isempty(rejected_comps_heart_PIERRE) == 0
            classification_iter(rejected_comps_heart_PIERRE) = 1;
        end
        % Gather will other recordings
        classification_heart_for_ROC = [classification_heart_for_ROC, classification_iter];



        % Gather all metrics values with other recordings
        % % % % % % % % % % % % % % % % % 
        FINAL_ALL_FEATURES_table{recording} = FINAL_ALL_FEATURES_table_iter;



        % Extract the recording info of the IC (will be useful to manually correct the classification for false negative and fasle positive)
        % % % % % % % % % % % % % % % % 
        IC_ori_iter = [];
        for IC_iter = 1:length(comp.label)
            IC_ori_iter = [IC_ori_iter; {[file_info '_IC_' num2str(IC_iter)]}];
        end
        IC_ori = [IC_ori; IC_ori_iter];


        % Gather some info( corr_matrix, z_scores) with the one of other recordings
        % % % % % % % % % % % % % % % % % 
        matrice_all{length(matrice_all)+1} = corr_matrice_iter;

        if length(rejected_comps_heart_PIERRE) == 1 % In case several heart IC are found
            zscore_ALL_heart_IC{1} = [zscore_ALL_heart_IC{1}; z_score_heart_IC.Properties.VariableNames];
            zscore_ALL_heart_IC{2} = [zscore_ALL_heart_IC{2}; table2array(z_score_heart_IC)];

            zscore_ALL_heart_IC_sort{1} = [zscore_ALL_heart_IC_sort{1}; z_score_heart_IC_sort.Properties.VariableNames];
            zscore_ALL_heart_IC_sort{2} = [zscore_ALL_heart_IC_sort{2}; table2array(z_score_heart_IC_sort)];           
        end


        if isempty(table_IC_cond_details) == 0
            bpm_all_IC = table_IC_cond_details.IC_bpm_all{1};
        end

            method_reject_cardiac_IC_parameters.Pierre_fct.method_reject_cardiac_IC_Pierre = cfg.method_chosen;
            method_reject_cardiac_IC_parameters.Pierre_fct.method_chosen = method_chosen;



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % B - METHOD EEGLAB (based on IC topography only)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


            % 0) For some subject (like sujet = 10), eeg_autocorr_fftw gives an error because the function 'resample' is not in Matlab path... (even if addpath is already at the beginning of the script...)

            % % % % % % % % % % % % % % % % % % % % % % % % %
            % addpath(strcat(main_path_PC, '/Toolbox/fieldtrip-20230118/fieldtrip-20230118/external/signal'));
            % % % % % % % % % % % % % % % % % % % % % % % % %


       



          


            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % D- PLOT FOR CHECKING
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            if strcmp(save_png, 'y')
                % 0) Extract IC timecourse
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                IC_timecourse = [];
                for i = 1:size(comp.trial,2)
                    IC_timecourse = [IC_timecourse, comp.trial{1,i}];
                end

                % 1) Plot heart IC (with ICLABEL)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isempty(rejected_comps_heart_EEGLAB) == 0
                    figure
                    IC_to_plot = rejected_comps_heart_EEGLAB;
                    plot([1:size(IC_timecourse,2)], IC_timecourse(IC_to_plot,:))
                    title(['Cardiac component iclabel (IC ' num2str(IC_to_plot) ')'])
                    %%%%%%%%%%% SAVE fig as .png %%%%%%%%%%%
                    filename = strcat(dirout, filesep, 'icalabel/', file_info, '_Rejected_IC_TimeCourse.png');
                    saveas(gcf, filename);
                    close(gcf);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                % 2) Plot heart IC (with PIERRE fct)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isempty(rejected_comps_heart_PIERRE) == 0
                    figure
                    IC_to_plot = rejected_comps_heart_PIERRE;
                    plot([1:size(IC_timecourse,2)], IC_timecourse(IC_to_plot,:))
                    title(['Cardiac component Pierre fct (IC ' num2str(IC_to_plot) ')'])
                    %%%%%%%%%%% SAVE fig as .png %%%%%%%%%%%
                    filename = strcat(dirout, filesep, 'Pierre_fct/', file_info, '_Rejected_IC_TimeCourse.png');
                    saveas(gcf, filename);
                    close(gcf);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end



            % %% V - 5) Identify eye IC (using SASICA)
            %
            % % 0) Detect bad cheek channels
            % % % % % % % % % % % % % % % % %
            %
            % % Keep cheek EEG
            % cfg = [];
            % cfg.channel = layout.label(cheeks_idx);
            % data_for_bad_cheek_EEG_detection = ft_preprocessing(cfg, data_cont_all_channels);
            %
            % % Identify bad channels within cheek EEG
            % data_for_bad_channel_detection = data_for_bad_cheek_EEG_detection;
            % run z_bad_channels_identification_PC.m;
            % bad_channels_name_cheek = bad_channels_name;
            %
            % bad_channels_for_EOG =[bad_channels_name_scalp; bad_channels_name_cheek];
            %
            % % 1) Define EEG channels that will be use to reconstruct EOG
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % V_EOG = {'E37', 'E241'};
            % H_EOG = {'E234', 'E244'};
            % bad_EEG_used_for_EOG = 0;
            %
            % % if strcmp(data_type, 'no_interp')
            % % (V_EOG) Try to find a combination of relevant EEG channels which have not been labelled as bad
            % if sum(ismember(V_EOG, bad_channels_for_EOG)) > 0
            %     V_EOG = {'E18', 'E238'};
            %     if sum(ismember(V_EOG, bad_channels_for_EOG)) > 0
            %         V_EOG = {'E46', 'E244'};
            %         if sum(ismember(V_EOG, bad_channels_for_EOG)) > 0
            %             V_EOG = {'E10', 'E234'};
            %             if sum(ismember(V_EOG, bad_channels_for_EOG)) > 0
            %                 V_EOG = {'E37', 'E241'};
            %                 bad_EEG_used_for_EOG = 1;
            %             end
            %         end
            %     end
            % end
            %
            % % (H_EOG) Try to find a combination of relevant EEG channels which have not been labelled as bad
            % if sum(ismember(H_EOG, bad_channels_for_EOG)) > 0
            %     H_EOG = {'E234', 'E244'};
            %     if sum(ismember(H_EOG, bad_channels_for_EOG)) > 0
            %         H_EOG = {'E239', 'E242'};
            %         if sum(ismember(H_EOG, bad_channels_for_EOG)) > 0
            %             H_EOG = {'E240', 'E243'};
            %             if sum(ismember(H_EOG, bad_channels_for_EOG)) > 0
            %                 H_EOG = {'E235', 'E245'};
            %                 if sum(ismember(H_EOG, bad_channels_for_EOG)) > 0
            %                     H_EOG = {'E234', 'E244'};
            %                     bad_EEG_used_for_EOG = 1;
            %                 end
            %             end
            %         end
            %     end
            % end
            %
            % % 3) Identify IC eye using SASICA (by correlating time course IC to EOG)
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % % addpath(strcat(path_toolbox, 'SASICA-master_PC/'));
            % cfg = [];
            % cfg.layout = layout;
            % cfg.EOGcorr.enable = 1;
            % cfg.EOGcorr.Veogchannames = V_EOG;
            % cfg.EOGcorr.Heogchannames = H_EOG;
            % cfg = ft_SASICA(cfg, comp, data_for_SASICA); %  / !! \ Use EEG dataset with all the 256 channels because we use the layout with the 256 channels
            %
            % % 4) Store the r-values of the correlation with EOG
            % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % sasica_output = cfg;
            % corr_vOEG = sasica_output.reject.SASICA.icachancorrVEOG;
            % corr_hEOG = sasica_output.reject.SASICA.icachancorrHEOG;
            %
            %
            % % 5) Identify the rejected IC
            % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % threshold_vEOG = mean(corr_vOEG) + IC_std_vEOG*std(corr_vOEG); %sasica_output.reject.SASICA.icathreshchancorrVEOG;
            % threshold_hEOG = mean(corr_hEOG) + IC_std_hEOG*std(corr_hEOG); %sasica_output.reject.SASICA.icathreshchancorrHEOG;
            %
            % rejected_comps_eye = unique([find(corr_vOEG>=threshold_vEOG & corr_vOEG>IC_r_min), find(corr_hEOG>=threshold_hEOG & corr_hEOG>IC_r_min)] ); %[find(sasica_output.reject.gcompreject == 1)];
            %
            %
            % % 6) Re-do the plot of r with EOG done by ft_SASICA and save it as .png
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %
            % [hplotcorr] = plot([corr_vOEG; corr_hEOG]','.','linestyle','none', 'MarkerSize', 15);
            % hold all
            % xlim([0 length(corr_hEOG)+1]);
            % xl = xlim;
            % yl = ylim;
            % set(gca,'ylimMode','manual');
            % plot(xl, [threshold_vEOG, threshold_vEOG], 'c--', 'LineWidth', 1);
            % plot(xl, [threshold_hEOG, threshold_hEOG], 'r--', 'LineWidth', 1);
            % title(['Correlation with EOG (vEOG: +' num2str(IC_std_vEOG) 'std, hEOG: + ' num2str(IC_std_hEOG) 'std, rmin:', num2str(IC_r_min) ')'])
            % ylabel('Correlation coef (r)');
            % xlabel('Components');
            %
            % %%%%% SAVE fig as .png %%%%%%%%
            % filename = strcat(dirout, filesep, file_info, '_EOG_IC_label (vEOG+', num2str(IC_std_vEOG), 'std,_hEOG+', num2str(IC_std_hEOG), 'std,_rmin', num2str(IC_r_min), ').png');
            % saveas(gcf, filename);
            % close(gcf);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %
            % % Gather all rejected components
            % rejected_comps = unique([rejected_comps_eye, rejected_comps_heart]);
            %
            %
            % % Plot rejected IC
            % if strcmp(save_png, 'y')
            %     if isempty(rejected_comps) == 0
            %         % 1) Look at IC topography + timecourse
            %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         cfg          = [];
            %         cfg.channel  = [rejected_comps]; % components to be plotted
            %         cfg.viewmode = 'component';
            %         cfg.layout   = layout_short;
            %
            %         cfg.trl(:, 1) = data_before_ICA.sampleinfo(epoch_to_display, 1);  % Start of each trial
            %         cfg.trl(:, 2) = data_before_ICA.sampleinfo(epoch_to_display, 2);  % End of each trial
            %         cfg.trl(:, 3) = 0;  % Offset
            %
            %         % If only 1 epoch (for FCSRT and WordList), plot only the first 30s
            %         if (cfg.trl(:, 2) - cfg.trl(:, 1) ) / data_before_ICA.fsample > 30
            %             cfg.trl(:, 2) = cfg.trl(:, 1) + 30*data_before_ICA.fsample;
            %         end
            %
            %         % Plot
            %         ft_databrowser(cfg, comp)
            %
            %         %%%%%%%%%%% SAVE fig as .png %%%%%%%%%%%
            %         filename = strcat(dirout, filesep, file_info, '_Rejected_IC_TimeCourse.png');
            %         saveas(gcf, filename);
            %         close(gcf);
            %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     end
            % end



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% VIII - Check visually ICA quality
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if strcmp(save_png, 'y')

                % 1) Plot data before ICA
                %%%%%%%%%%%%%%%%%%%%%%%%%
                cfg = [];
                cfg.continuous  = 'no';
                cfg.viewmode    = 'vertical'; % all channels separate
                cfg.blocksize   = 30;         % view the continuous data in 30-s epochs
                % cfg.artfctdef.art.artifact = art_FINAL_trl;
                cfg.channel = data_before_ICA.label([1:6, 25, 26, 31, 40, 150:160]);
                % Display epoch 2 (if RS) because subjects have eyes closed in epoch 1 during RS so no blinks (to compare with IC time course)
                cfg.trl(:, 1) = data_before_ICA.sampleinfo(epoch_to_display, 1);  % Start of each trial
                cfg.trl(:, 2) = data_before_ICA.sampleinfo(epoch_to_display, 2);  % End of each trial
                cfg.trl(:, 3) = 0;  % Offset

                ft_databrowser(cfg, data_before_ICA);

                %%%%%%%%%%% SAVE fig as .png %%%%%%%%%%%
                filename = strcat(dirout, filesep, file_info, '_EEG_Epoch', num2str(epoch_to_display), '_Before_ICA.png');
                saveas(gcf, filename);
                close(gcf);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                % 2) Plot data after ICA
                %%%%%%%%%%%%%%%%%%%%%
                cfg = [];
                cfg.continuous  = 'no';
                cfg.viewmode    = 'vertical'; % all channels separate
                cfg.blocksize   = 30;         % view the continuous data in 30-s epochs
                % cfg.artfctdef.art.artifact = art_FINAL_trl;
                cfg.channel = data_clean_not_filt.label([1:6, 25, 26, 31, 40, 150:160]);
                % Display epoch 2 (if RS) because subjects have eyes closed in epoch 1 during RS so no blinks (to compare with IC time course)
                cfg.trl(:, 1) = data_clean_not_filt.sampleinfo(epoch_to_display, 1);  % Start of each trial
                cfg.trl(:, 2) = data_clean_not_filt.sampleinfo(epoch_to_display, 2);  % End of each trial
                cfg.trl(:, 3) = 0;  % Offset

                ft_databrowser(cfg, data_clean_not_filt);

                %%%%%%%%%%% SAVE fig as .png %%%%%%%%%%%
                filename = strcat(dirout, filesep, file_info, '_EEG_Epoch', num2str(epoch_to_display), '_After_ICA.png');
                saveas(gcf, filename);
                close(gcf);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end






            %% XI - Restore path (avoid issues using FieldTrip + SASICA that both contain different versions and functions of EEGLAB)

            restoredefaultpath
            run z_paths_LS_Insight_PC.m

        end
    end
end


%% ROC analysis


% 0) USER manually correct classification_heart_for_ROC (use IC_ori)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_false_neg = {'S004_M0_run_1_IC_10', 'S004_M0_run_7_IC_12', ...
'S023_M60_run_1_IC_2',  'S023_M60_run_7_IC_2', 'S060_M24_run_1_IC_13', 'S098_M36_run_1_IC_3', ...
'S100_M60_run_1_IC_10', 'S203_M24_run_1_IC_2', 'S203_M24_run_7_IC_16', 'S314_M0_run_1_IC_18', ...
'S314_M12_run_1_IC_21', 'S314_M48_run_1_IC_7', 'S314_M48_run_7_IC_5',  'S314_M60_run_7_IC_8', ...
'S368_M0_run_1_IC_1',   'S368_M0_run_7_IC_1',  'S368_M36_run_1_IC_2',  'S368_M36_run_7_IC_3', ...
'S368_M60_run_1_IC_4',  'S370_M24_run_1_IC_4', 'S370_M36_run_1_IC_5',  'S370_M36_run_7_IC_5'};

idx_false_neg = find(ismember(IC_ori, list_false_neg));

classification_heart_for_ROC(idx_false_neg) = 1;

% 1) Prepare values of the metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USER: Chose metric for ROC
var_of_interest = 'ERPampl_median';

% Find column of var interest
var_list = FINAL_ALL_FEATURES_table_iter.Properties.VariableNames;
var_idx = find(ismember(var_list, var_of_interest));

% Extract the values for all the IC of all the recordings (for var of interest only)
score_var_interest = [];
for recording_iter = 1:length(FINAL_ALL_FEATURES_table)
    data = FINAL_ALL_FEATURES_table{recording_iter};

    score_var_interest_iter = data(:,var_idx);

    score_var_interest = [score_var_interest, table2array(score_var_interest_iter)'];
end

% Change score sign bc we must have cardiac predicted by high value for ROC
score_var_interest = - score_var_interest; 


% 2) Prepare the true classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classification_heart_for_ROC = classification_heart_for_ROC';



% 3) ROC Curve
%%%%%%%%%%%%
[X, Y, T, AUC] = perfcurve(classification_heart_for_ROC', score_var_interest', 1);

% Plot
plot(X, Y)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title(['ROC Curve (AUC = ' num2str(AUC) ')'])







%% ALL ROC CURVES ON SAME PLOT


% 0) USER manually correct classification_heart_for_ROC (use IC_ori)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_false_neg = {'S004_M0_run_1_IC_10', 'S004_M0_run_7_IC_12', ...
'S023_M60_run_1_IC_2',  'S023_M60_run_7_IC_2', 'S060_M24_run_1_IC_13', 'S098_M36_run_1_IC_3', ...
'S100_M60_run_1_IC_10', 'S203_M24_run_1_IC_2', 'S203_M24_run_7_IC_16', 'S314_M0_run_1_IC_18', ...
'S314_M12_run_1_IC_21', 'S314_M48_run_1_IC_7', 'S314_M48_run_7_IC_5',  'S314_M60_run_7_IC_8', ...
'S368_M0_run_1_IC_1',   'S368_M0_run_7_IC_1',  'S368_M36_run_1_IC_2',  'S368_M36_run_7_IC_3', ...
'S368_M60_run_1_IC_4',  'S370_M24_run_1_IC_4', 'S370_M36_run_1_IC_5',  'S370_M36_run_7_IC_5'};

idx_false_neg = find(ismember(IC_ori, list_false_neg));

classification_heart_for_ROC(idx_false_neg) = 1;



var_list = FINAL_ALL_FEATURES_table_iter.Properties.VariableNames;

classification_heart_for_ROC = classification_heart_for_ROC';

% Option 1: All ROC curves on same figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AUC_all = [];
X_MRMR = [];

figure;
for i = 1:length(var_list)
    var_idx = i;

    % Extract the values for all the IC of all the recordings (for var of interest only)
    score_var_interest = [];
    for recording_iter = 1:length(FINAL_ALL_FEATURES_table)
        data = FINAL_ALL_FEATURES_table{recording_iter};

        score_var_interest_iter = data(:,var_idx);

        score_var_interest = [score_var_interest, table2array(score_var_interest_iter)'];
    end

    % Change score sign bc we must have cardiac predicted by high value for ROC
    score_var_interest = - score_var_interest;

    % Gather with other metrics for MRMR
    X_MRMR = [X_MRMR, score_var_interest']


    % ROC Computation
    [X, Y, T, AUC] = perfcurve(classification_heart_for_ROC, score_var_interest, 1);

    % Gather AUC with other metrics
    AUC_all = [AUC_all; data.Properties.VariableNames(var_idx), AUC];

    % Plot
    hold on;
    plot(X, Y)
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    % title(['ROC Curve (AUC = ' num2str(AUC) ')'])
    hold on;

end

% Plot AUC by decreasing values
AUC_all_sorted = sortrows(AUC_all, 2, 'descend');

names = AUC_all_sorted(:,1);
values = cell2mat(AUC_all_sorted(:,2));

figure;
bar(values); 
xticks(1:length(names));           % Set the tick positions
xticklabels(names);                % Set the tick labels to the names
xtickangle(45);                    % Rotate x-axis labels for better readability
ylabel('AUC');

% Option 2: A figure per ROC curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER: Chose metric for ROC
var_of_interest = 'Sampl_ratio_std_mean';
var_idx = find(ismember(var_list, var_of_interest));

% Extract the values for all the IC of all the recordings (for var of interest only)
score_var_interest = [];
for recording_iter = 1:length(FINAL_ALL_FEATURES_table)
    data = FINAL_ALL_FEATURES_table{recording_iter};

    score_var_interest_iter = data(:,var_idx);

    score_var_interest = [score_var_interest, table2array(score_var_interest_iter)'];
end

% Change score sign bc we must have cardiac predicted by high value for ROC
score_var_interest = - score_var_interest; 


% ROC computation
[X, Y, T, AUC] = perfcurve(classification_heart_for_ROC', score_var_interest', 1);

% Plot
figure;
plot(X, Y)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title(['ROC Curve (AUC = ' num2str(AUC) ')'])


%% MRMR

idx_var = 1:9


% Extract feature names
list_metric = FINAL_ALL_FEATURES_table{1, 1}  ;
list_metric = list_metric.Properties.VariableNames(idx_var);

% Prepare MRMR
Y_MRMR = classification_heart_for_ROC';

% Run MRMR
[feature_rank, feature_score] = fscmrmr(X_MRMR(:,idx_var), Y_MRMR);

MRMR_output = [list_metric(idx_var)', num2cell(feature_score')];

% Sort features by decreasing MRMR score
% MRMR_output = sortrows(MRMR_output, 2, 'descend');

% Plot
names = MRMR_output(:,1);
values = cell2mat(MRMR_output(:,2));

figure;
bar(values); 
xticks(1:length(names));           % Set the tick positions
xticklabels(names);                % Set the tick labels to the names
xtickangle(45);                    % Rotate x-axis labels for better readability
ylabel('MRMR score');

%% MRMR with balanced dataset

idx_var = 1:9


% Extract feature names
list_metric = FINAL_ALL_FEATURES_table{1, 1}  ;
list_metric = list_metric.Properties.VariableNames(idx_var);

% Prepare MRMR
Y_MRMR = classification_heart_for_ROC';

% Randomly eliminate non-cardiac IC

idx_to_keep = randperm(length(find(Y_MRMR==0)), length(find(Y_MRMR==1)));
idx_to_keep = [idx_to_keep, find(Y_MRMR==1)]

X_MRMR_separate = X_MRMR(idx_to_keep,idx_var);
Y_MRMR_separate = Y_MRMR(idx_to_keep);

% Run MRMR
[feature_rank, feature_score] = fscmrmr(X_MRMR_separate, Y_MRMR_separate);

MRMR_output = [list_metric', num2cell(feature_score')];

% Sort features by decreasing MRMR score
% MRMR_output = sortrows(MRMR_output, 2, 'descend');

% Plot
names = MRMR_output(:,1);
values = cell2mat(MRMR_output(:,2));

figure;
bar(values); 
xticks(1:length(names));           % Set the tick positions
xticklabels(names);                % Set the tick labels to the names
xtickangle(45);                    % Rotate x-axis labels for better readability
ylabel('MRMR score');



% combined AUC 
y = Y_MRMR_separate;
AUC_for_loop = [];

for i = 1:size(feature_rank,2)
    X = X_MRMR_separate(:, feature_rank(1:i));


    % Combine metrics using logistic regression
    mdl = fitglm(X, y, 'Distribution', 'binomial', 'Link', 'logit');

    % Predict combined scores
    scores = predict(mdl, X);

    % Compute ROC curve
    [X_roc, Y_roc, T, AUC] = perfcurve(y, scores, 1);

    % Gather
    AUC_for_loop_iter = [{list_metric(feature_rank(1:i))}, AUC];
    AUC_for_loop = [AUC_for_loop; AUC_for_loop_iter];

end

figure;
plot(1:length(idx_var), cell2mat(AUC_for_loop(:,2)));
ylabel('combined AUC (logistic regression)')
xlabel('Nb of features combined')


%% AUC of combined metrics

X = X_MRMR(:, feature_rank(1:3));
y = classification_heart_for_ROC;

% Combine metrics using logistic regression
mdl = fitglm(X, y, 'Distribution', 'binomial', 'Link', 'logit');

% Predict combined scores
scores = predict(mdl, X);

% Compute ROC curve
[X_roc, Y_roc, T, AUC] = perfcurve(y, scores, 1);

% Plot the ROC curve
figure;
plot(X_roc, Y_roc, 'LineWidth', 2);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(['ROC Curve (AUC = ' num2str(AUC) ')']);
grid on;



%%%%%%%%%%%%%%%%%%%%%%%% WITH FOR LOOP

y = classification_heart_for_ROC;
AUC_for_loop = [];

for i = 1:size(feature_rank,2)
    X = X_MRMR(:, feature_rank(1:i));


    % Combine metrics using logistic regression
    mdl = fitglm(X, y, 'Distribution', 'binomial', 'Link', 'logit');

    % Predict combined scores
    scores = predict(mdl, X);

    % Compute ROC curve
    [X_roc, Y_roc, T, AUC] = perfcurve(y, scores, 1);

    % Gather
    AUC_for_loop_iter = [{list_metric(feature_rank(1:i))}, AUC];
    AUC_for_loop = [AUC_for_loop; AUC_for_loop_iter];

end

figure;
plot([1:length(feature_rank)], cell2mat(AUC_for_loop(:,2)));
ylabel('combined AUC (logistic regression)')
xlabel('Nb of features combined')


%% Random Forest

% Inputs 
idx_var = 1:9
X = X_MRMR(:,idx_var);
Y = classification_heart_for_ROC;


% Calculate class weights
classWeights = ones(size(Y));
classWeights(Y == 1) = length(Y) / sum(Y == 1);
classWeights(Y == 0) = length(Y) / sum(Y == 0);


% Train Random Forest with class weights
RF_model = TreeBagger(100, X, Y, 'Method', 'classification', 'Weights', classWeights, 'OOBPredictorImportance', 'on');


% Get feature importance
importance = RF_model.OOBPermutedVarDeltaError;
disp('Feature Importance:');
disp(importance);


[~, sorted_indices] = sort(importance, 'descend');  % Sort importance
top_features = sorted_indices(1:size(X,2));  % Get indices of top 5 features
list_metric(top_features)


RF_output = [list_metric(idx_var)', num2cell(importance)'];

% Plot features by decreasing MRMR score
% RF_output = sortrows(RF_output, 2, 'descend');

names = RF_output(:,1);
values = cell2mat(RF_output(:,2));

figure;
bar(values); 
xticks(1:length(names));           % Set the tick positions
xticklabels(names);                % Set the tick labels to the names
xtickangle(45);                    % Rotate x-axis labels for better readability
ylabel('RF score');


%% Plot and save an average corr matrix

meanMatrix = mean(cat(3, matrice_all{:}), 3); % Concatenate along the 3rd dimension and take mean

min_scale = -1;
max_scale = 1;

% List variable names for corr matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varNames = {'ERPampl_median', 'ERPampl_std', 'ERPampl_skew', 'ERPampl_kurt', 'ERPampl_std_median',...
    'RR_ratio_std_mean', 'RR_skew', 'RR_kurt', 'RR_chi2stat', ...
    'Rampl_ratio_std_mean', 'Rampl_skew', 'Rampl_kurt', 'Rampl_chi2stat', ...
    'Qampl_ratio_std_mean', 'Qampl_skew', 'Qampl_kurt', 'Qampl_chi2stat', ...
    'Sampl_ratio_std_mean', 'Sampl_skew', 'Sampl_kurt', 'Sampl_chi2stat'};


% Plot corr matrix
%%%%%%%%%%%%%%%%%%
figure
imagesc(meanMatrix);
colorbar;
colormap(jet);
xticks(1:length(varNames)); xticklabels(varNames);
yticks(1:length(varNames)); yticklabels(varNames);
% Set color limits (e.g., [-1, 1] for correlation)
caxis([min_scale max_scale]);                          % Define color bar limits
colorbar;                               % Add color bar
colormap(jet);                          % Use a jet colormap
title(['Correlation Matrix (AVERAGE ' num2str(length(matrice_all)) ' RECORDINGS)']);


% Annotate each cell with its r-value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [nRows, nCols] = size(correlationMatrix); % Get the matrix dimensions
% for i = 1:nRows
%     for j = 1:nCols
%         text(j, i, num2str(correlationMatrix(i, j), '%.2f'), ... % Format values to 2 decimal places
%             'HorizontalAlignment', 'center', ...
%             'Color', 'white', ...         % Text color
%             'FontSize', 10);             % Adjust font size
%     end
% end

% Save
%%%%%%
% filename = strcat([cfg.path_output, '/corr_matrices/AA_AVERAGE ' num2str(length(matrice_all)) ' RECORDINGS_Corr_matrix.png']);
% saveas(gcf, filename);
% close(gcf);



%% Which feature is most useful for zscore

best_zscore_nb = '3';
keep_ERP_ampl = 'y'; % y/n

% Exclude column of ERP_ampl if needed
if strcmp(keep_ERP_ampl, 'n')
    zscore_ALL_heart_IC_sort_NEW = [];
    for i = 1:size(zscore_ALL_heart_IC_sort{1, 1}, 1)
        idx_ERP = find(ismember(zscore_ALL_heart_IC_sort{1, 1}(i,:), {'z_ERPampl'}));
        zscore_ALL_heart_IC_sort_NEW = [zscore_ALL_heart_IC_sort_NEW; zscore_ALL_heart_IC_sort{1, 1}(i, [1:idx_ERP-1, idx_ERP+1:end])];
    end

    data = zscore_ALL_heart_IC_sort_NEW;
else
    data = zscore_ALL_heart_IC_sort{1, 1};
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Zscore n°1 to n° "best_zscore_nb"
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if strcmp(best_zscore_nb, '1') == 0 
     % Select data column 1 to column "best_zscore_nb"
    data_to_plot = categorical(data(:,1:str2num(best_zscore_nb)));
    data_to_plot = addcats(data_to_plot, data(1,:));

    % Reorder categorical levels
    if strcmp(keep_ERP_ampl, 'n')
        data_to_plot = reordercats(data_to_plot, {'z_RR_ratio_std_mean', 'z_RR_skew', 'z_RR_kurt', 'z_RR_chi2stat', ...
            'z_Rampl_ratio_std_mean', 'z_Rampl_skew', 'z_Rampl_kurt', 'z_Rampl_chi2stat', ...
            'z_Qampl_ratio_std_mean', 'z_Qampl_skew', 'z_Qampl_kurt', 'z_Qampl_chi2stat', ...
            'z_Sampl_ratio_std_mean', 'z_Sampl_skew', 'z_Sampl_kurt', 'z_Sampl_chi2stat'});
    else
        data_to_plot = reordercats(data_to_plot, {'z_ERPampl_median', 'z_ERPampl_std', 'z_ERPampl_skew', 'z_ERPampl_kurt', 'z_ERPampl_std_median', 'z_RR_ratio_std_mean', 'z_RR_skew', 'z_RR_kurt', 'z_RR_chi2stat', ...
            'z_Rampl_ratio_std_mean', 'z_Rampl_skew', 'z_Rampl_kurt', 'z_Rampl_chi2stat', ...
            'z_Qampl_ratio_std_mean', 'z_Qampl_skew', 'z_Qampl_kurt', 'z_Qampl_chi2stat', ...
            'z_Sampl_ratio_std_mean', 'z_Sampl_skew', 'z_Sampl_kurt', 'z_Sampl_chi2stat'});
    end
    
    % Plot
    figure
    histogram(data_to_plot, 'Normalization', 'probability')
    yticks = get(gca, 'YTick');
    set(gca, 'YTickLabel', yticks * 100);
    ylabel('% of recordings')
    title(['Best zscore n° 1 to ' best_zscore_nb ' to identify heart IC (' num2str(size(data,1)) ' recordings)'])
 
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Only Zscore n° "best_zscore_nb"
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Select data column 1 to column "best_zscore_nb"
data_to_plot = categorical(data(:,str2num(best_zscore_nb)));
data_to_plot = addcats(data_to_plot, data(1,:));

% Reorder categorical levels
if strcmp(keep_ERP_ampl, 'n')
    data_to_plot = reordercats(data_to_plot, {'z_RR_ratio_std_mean', 'z_RR_skew', 'z_RR_kurt', 'z_RR_chi2stat', ...
            'z_Rampl_ratio_std_mean', 'z_Rampl_skew', 'z_Rampl_kurt', 'z_Rampl_chi2stat', ...
            'z_Qampl_ratio_std_mean', 'z_Qampl_skew', 'z_Qampl_kurt', 'z_Qampl_chi2stat', ...
            'z_Sampl_ratio_std_mean', 'z_Sampl_skew', 'z_Sampl_kurt', 'z_Sampl_chi2stat'});
else
    data_to_plot = reordercats(data_to_plot, {'z_ERPampl_median', 'z_ERPampl_std', 'z_ERPampl_skew', 'z_ERPampl_kurt', 'z_ERPampl_std_median', 'z_RR_ratio_std_mean', 'z_RR_skew', 'z_RR_kurt', 'z_RR_chi2stat', ...
            'z_Rampl_ratio_std_mean', 'z_Rampl_skew', 'z_Rampl_kurt', 'z_Rampl_chi2stat', ...
            'z_Qampl_ratio_std_mean', 'z_Qampl_skew', 'z_Qampl_kurt', 'z_Qampl_chi2stat', ...
            'z_Sampl_ratio_std_mean', 'z_Sampl_skew', 'z_Sampl_kurt', 'z_Sampl_chi2stat'});
end

% Plot
figure
histogram(data_to_plot, 'Normalization', 'probability')
yticks = get(gca, 'YTick');        
set(gca, 'YTickLabel', yticks * 100);
ylabel('% of recordings')
title(['Best zscore n°' best_zscore_nb ' to identify heart IC (' num2str(size(data,1)) ' recordings)'])

% save([pwd '/Matlab_workspace_zscore_12conditions.mat'])


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Quantification of zscore values for each metric
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

figure
bar(categorical(zscore_ALL_heart_IC{1,1}(1,:)), median(zscore_ALL_heart_IC{1,2}))

% yticks = get(gca, 'YTick');        
% set(gca, 'YTickLabel', yticks * 100);
ylabel('Median z score')
title(['Median zscore in cardiac IC (' num2str(size(data,1)) ' recordings)'])
