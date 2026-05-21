
which_data = '_SASICA_replicate'; % '' // '_30comp'_HPF1Hz

setpath_ds003690

compare_truth = 'MANUAL';
compare_with = 'SASICARACAS'; %'ICLabel'; % 'CARACAS'; % 'CORR';

rm_frompath('eeglab')

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

rng(123)
jitt = 0.1;

outerK = 5;
CV = cvpartition(length(fs), 'KFold', outerK);
outer_results = zeros(outerK, 1); 

% We store each outer fold's table in a cell array to avoid BigT collision
BigT_cells = cell(outerK, 1);

for i_train = 1:outerK
    Trainfs = fs(CV.training(i_train));
    Testfs = fs(CV.test(i_train));
    
    %% Setup Initial Guess
    param_vec = [1, 5, 1/3, 1/4, 35, 100];
    param_vec = param_vec .* (1 + jitt * randn(size(param_vec)));
    
    %% Inner K-fold (Parallelized)
    innerK = 5;
    cv_inner = cvpartition(length(Trainfs), 'KFold', innerK);
    
    % Pre-allocate for parfor
    inner_best_params = zeros(innerK, length(param_vec));
    inner_errors = zeros(innerK, 1);
    innerSensitivity = zeros(innerK, 1);
    innerSpecificity = zeros(innerK, 1);
    
    % Use parfor for the optimization heavy-lifting
    parfor k = 1:innerK
        % Extract indices inside the loop
        trainIdx = training(cv_inner, k);
        testIdx  = test(cv_inner, k);
        
        % Local function handle for workers
        myoptim = @(X) CARACASerr(X, Trainfs(trainIdx), compare_truth, compare_with);
        
        % Run optimization
        inner_best_params(k,:) = fminsearch(myoptim, param_vec, struct('Display','off'));
        [inner_errors(k), innerSensitivity(k), innerSpecificity(k)] = CARACASerr(inner_best_params(k,:), Trainfs(testIdx), compare_truth, compare_with);
    end
    
    %% Prepare Table for this Outer Fold (Post-Parallel processing)
    param_names = {'Start_sk', 'Start_ku', 'Start_RR', 'Start_Rampl', 'Start_bpm1', 'Start_bpm2', ...
                   'Best_sk', 'Best_ku', 'Best_RR', 'Best_Rampl', 'Best_bpm1', 'Best_bpm2', 'Error', 'Sensitivity', 'Specificity'};
    
    % Prepare the rows for this fold
    all_fold_data = [repmat(param_vec, innerK, 1), inner_best_params, inner_errors, innerSensitivity, innerSpecificity];
    mean_params = mean(inner_best_params, 1);
    mean_err = mean(inner_errors);
    meanSens = mean(innerSensitivity);
    meanSpec = mean(innerSpecificity);
    
    % Final Test for this outer slice
    [finalTestErr, finalSensitivity, finalSpecificity] = CARACASerr(mean_params, Testfs, compare_truth, compare_with);
    outer_results(i_train) = finalTestErr;
    
    % Assemble data for current outer fold
    fold_matrix = [all_fold_data; 
                   param_vec, mean_params, mean_err, meanSens, meanSpec; 
                   param_vec, mean_params, finalTestErr, finalSensitivity, finalSpecificity];
               
    T = array2table(fold_matrix, 'VariableNames', param_names);
    
    % Set Row Names
    rowLabs = [arrayfun(@(x) sprintf('Out%d_In%d', i_train, x), 1:innerK, 'UniformOutput', false), ...
               {sprintf('Out%d_Mean', i_train)}, {sprintf('Out%d_FINAL', i_train)}];
    T.Properties.RowNames = rowLabs;
    
    BigT_cells{i_train} = T;
    
    fprintf('Finished Outer Fold %d/%d. Test Error: %.4g Sensitivity: %.4g Specificity: %.4g\n', i_train, outerK, finalTestErr, finalSensitivity, finalSpecificity);
end

% Concatenate all tables at the end
BigT = vertcat(BigT_cells{:});

%% Final Reporting
disp('--- Full Nested CV Results ---')
disp(BigT)
fprintf('Grand Mean Error: %.4f\n', mean(outer_results));
fprintf('Standard Deviation: %.4f\n', std(outer_results));

writetable(BigT,sprintf('script_05_optim_recross_%g.csv',jitt), WriteRowNames=true)

function [myerr, Sensitivity, Specificity] = CARACASerr(threshs, fs, compare_truth, compare_with)

cfg_SASICA = SASICA('getdefs');
cfg_CARACAS = cfg_SASICA.CARACAS;
cfg_CARACAS.thresh_sk = threshs(1);
cfg_CARACAS.thresh_ku = threshs(2);
cfg_CARACAS.thresh_RR = threshs(3);
cfg_CARACAS.thresh_Rampl = threshs(4);
cfg_CARACAS.thresh_bpm = threshs(5:6);

switch compare_truth
    case 'MANUAL'
        for i = 1:numel(fs)
            MANUALrej(i,:) = fs(i).MANUAL.rej;
        end
        truthrej = MANUALrej;
    case 'SASICARACAS'
        SASICARACASrej = zeros(numel(fs),numel(fs(1).SASICARARAS.rej));
        for i = 1:numel(fs)
            SASICARACASrej(i,:) = CARACAS_rethresh(SASICARACASrej(i,:), fs(i).SASICARACAS, cfg_CARACAS);
        end
        truthrej = SASICARACASrej;
    case 'CARACAS'
        CARACASrej = zeros(numel(fs),numel(fs(1).CARARAS.rej));
        for i = 1:numel(fs)
            CORRrej(i,:) = fs(i).CORR.rej;

            CARACASrej(i,:) = CARACAS_rethresh(CARACASrej(i,:), fs(i).CARACAS, cfg_CARACAS);

            SASICARACASrej(i,:) = CARACAS_rethresh(SASICARACASrej(i,:), fs(i).SASICARACAS, cfg_CARACAS);
            MANUALrej(i,:) = fs(i).MANUAL.rej;
        end
        truthrej = CARACASrej;
    case 'CORR'
        CORRrej = zeros(numel(fs),numel(fs(1).CORR.rej));
        for i = 1:numel(fs)
            CORRrej(i,:) = fs(i).CORR.rej;
        end
        truthrej = CORRrej;
end
switch compare_with
    case 'MANUAL'
        for i = 1:numel(fs)
            MANUALrej(i,:) = fs(i).MANUAL.rej;
        end
        withrej = MANUALrej;
    case 'SASICARACAS'
        SASICARACASrej = zeros(numel(fs),numel(fs(1).SASICARACAS.rej));
        for i = 1:numel(fs)
            SASICARACASrej(i,:) = CARACAS_rethresh(SASICARACASrej(i,:), fs(i).SASICARACAS, cfg_CARACAS);
        end
        withrej = SASICARACASrej;
    case 'CARACAS'
        CARACASrej = zeros(numel(fs),numel(fs(1).CARARAS.rej));
        for i = 1:numel(fs)
            CARACASrej(i,:) = CARACAS_rethresh(CARACASrej(i,:), fs(i).CARACAS, cfg_CARACAS);
        end
        withrej = CARACASrej;
    case 'CORR'
        CORRrej = zeros(numel(fs),numel(fs(1).CORR.rej));
        for i = 1:numel(fs)
            CORRrej(i,:) = fs(i).CORR.rej;
        end
        withrej = CORRrej;
end

toplot = NaN(size(truthrej));
toplot(truthrej & withrej) = 1;     % Hit
toplot(truthrej & ~ withrej) = 2;   % Miss
toplot(~truthrej & ~withrej) = 3;   % CR
toplot(~truthrej & withrej) = 4;    % FA
names = {'Hit','Miss','CR','FA'};

% Performances
Sensitivity = sum(toplot == 1) / (sum(toplot== 1) + sum(toplot == 2));
Specificity = sum(toplot == 3) / (sum(toplot== 3) + sum(toplot == 4));


Balanced_accuracy = (Sensitivity + Specificity) / 2;
% Weighted_accuracy = 0.25*Sensitivity + 0.75*Specificity;
% Accuracy = ( sum(toplot == 1) + sum(toplot == 3) ) / (sum(toplot == 1) + sum(toplot == 2) + sum(toplot == 3) + sum(toplot == 4) );

% myerr = 1 - Weighted_accuracy;
myerr = 1 - Balanced_accuracy;

end
