%% --- 1. SETUP & DATA LOADING ---
which_data = '_SASICA_replicate'; 
setpath_ds003690
compare_truth = 'MANUAL';
compare_with = 'SASICARACAS'; 
rm_frompath('eeglab')
load(fullfile(dirout,'AllFilesAndScoresList.mat'))

% Fix paths for PC if necessary
if ispc
    for i = 1:numel(fs)
        fields = {'name', 'compf','eegf'};
        for j = 1:numel(fields)
            fs(i).(fields{j}) = strrep(fs(i).(fields{j}),'/network','\');
            fs(i).(fields{j}) = strrep(fs(i).(fields{j}),'/','\');
        end
    end
end

%% --- 2. BENCHMARK CONFIGURATION ---
% Define the "Golden" starting point
base_param_vec = [1, 5, 1/3, 1/4, 35, 100];
param_labels = {'sk', 'ku', 'RR', 'Rampl', 'bpm1', 'bpm2'};

% Define the grid: 5 values from -20% to +20% (0.8x to 1.2x)
grid_multipliers = linspace(0.8, 1.2, 5); 

% Storage for global summary
% Columns: Param_Idx, Param_Name, Multiplier, Initial_Val, Mean_Err, Std_Err
benchmark_results = table();

fprintf('Starting Sensitivity Analysis: 6 parameters x 5 initial values...\n');

%% --- 3. BENCHMARK LOOP ---
for p_idx = 1:length(base_param_vec)
    for g_idx = 1:length(grid_multipliers)
        
        % Current multiplier
        mult = grid_multipliers(g_idx);
        
        % Prepare current initial vector
        param_vec = base_param_vec;
        param_vec(p_idx) = base_param_vec(p_idx) * mult;
        
        fprintf('Testing %s with initial value %.4f (%.1f%% of base)\n', ...
            param_labels{p_idx}, param_vec(p_idx), mult*100);

        % --- NESTED CV BLOCK (Your original logic) ---
        rng(123) % Fixed seed for reproducibility across benchmarks
        outerK = 5;
        CV = cvpartition(length(fs), 'KFold', outerK);
        outer_results = zeros(outerK, 1); 

        for i_train = 1:outerK
            Trainfs = fs(CV.training(i_train));
            Testfs = fs(CV.test(i_train));
            
            % Inner K-fold
            innerK = 5;
            cv_inner = cvpartition(length(Trainfs), 'KFold', innerK);
            inner_best_params = zeros(innerK, length(param_vec));
            inner_errors = zeros(innerK, 1);
            
            parfor k = 1:innerK
                trainIdx = training(cv_inner, k);
                testIdx  = test(cv_inner, k);
                myoptim = @(X) CARACASerr(X, Trainfs(trainIdx), compare_truth, compare_with);
                
                % Optimization
                inner_best_params(k,:) = fminsearch(myoptim, param_vec, struct('Display','off'));
                inner_errors(k) = CARACASerr(inner_best_params(k,:), Trainfs(testIdx), compare_truth, compare_with);
            end
            
            % Final Test for this outer slice
            mean_params = mean(inner_best_params, 1);
            outer_results(i_train) = CARACASerr(mean_params, Testfs, compare_truth, compare_with);
        end
        % --- END NESTED CV BLOCK ---

        % Append to summary table
        res_row = {p_idx, param_labels{p_idx}, mult, param_vec(p_idx), mean(outer_results), std(outer_results)};
        benchmark_results = [benchmark_results; res_row];
    end
end

% Set table variable names
benchmark_results.Properties.VariableNames = {'ParamIdx', 'ParamName', 'Multiplier', 'InitialVal', 'MeanErr', 'StdErr'};

%% --- 4. REPORTING ---
disp('--- Sensitivity Analysis Summary ---')
disp(benchmark_results)
writetable(benchmark_results, 'Sensitivity_Analysis_Results.csv')

% Optional: Quick plot to see trends
figure;
for p = 1:6
    subplot(2,3,p);
    idx = benchmark_results.ParamIdx == p;
    plot(benchmark_results.InitialVal(idx), benchmark_results.MeanErr(idx), '-o', 'LineWidth', 2);
    title(['Effect of ', param_labels{p}]);
    xlabel('Initial Start Value'); ylabel('Mean CV Error');
    grid on;
end



function myerr = CARACASerr(threshs, fs, compare_truth, compare_with)

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
