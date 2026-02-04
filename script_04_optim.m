
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

%% Number of folds
K = 5;
cfg_SASICA = SASICA('getdefs');
% cfg_SASICA.CARACAS.enable = 1;
% cfg_CARACAS = cfg_SASICA.CARACAS;
% cfg_CARACAS.thresh_RPeakstoNoise = 10;
% cfg_CARACAS.thresh_sk = 1;
% cfg_CARACAS.thresh_ku = 5;
% cfg_CARACAS.thresh_RR  = 1/3;
% cfg_CARACAS.thresh_Rampl = 1/4;
% cfg_CARACAS.thresh_bpm = [35 100];
param_vec = [cfg_CARACAS.thresh_RPeakstoNoise cfg_CARACAS.thresh_sk cfg_CARACAS.thresh_ku cfg_CARACAS.thresh_RR cfg_CARACAS.thresh_Rampl cfg_CARACAS.thresh_bpm];

% Split data into K folds
cv = cvpartition(length(fs), 'KFold', K);
% Store results
errors = zeros(K, 1);
params = zeros(K, 1);

testerr = [];bestparams = [];
for k = 1:K

   % Training and validation indices
    trainIdx = training(cv, k);
    testIdx  = test(cv, k);

    myoptim = @(X) CARACASerr(X,fs(trainIdx), compare_truth, compare_with);
    
    bestparams(k,:) = fminsearch(myoptim,param_vec, struct('Display','iter'));
% 
    testerr(k) = CARACASerr(bestparams(k,:),fs(testIdx), compare_truth, compare_with);
end



%%
% Display results as a table
param_names = {'thresh_RPeakstoNoise', 'thresh_sk', 'thresh_ku', 'thresh_RR', 'thresh_Rampl', 'thresh_bpm1', 'thresh_bpm2', 'TestError'};
results_table = [bestparams, testerr'];
mean_row = [mean(bestparams), mean(testerr)];

% Create table
T = array2table([results_table; mean_row], 'VariableNames', param_names);
T.Properties.RowNames = [arrayfun(@(x) sprintf('Fold_%d', x), 1:K, 'UniformOutput', false), {'Mean'}];

% Display table
disp('Cross-validation Results:')
disp(T)

%

function myerr = CARACASerr(threshs, fs, compare_truth, compare_with)

cfg_SASICA = SASICA('getdefs');
cfg_CARACAS = cfg_SASICA.CARACAS;
cfg_CARACAS.thresh_RPeakstoNoise = threshs(1);
cfg_CARACAS.thresh_sk = threshs(2);
cfg_CARACAS.thresh_ku = threshs(3);
cfg_CARACAS.thresh_RR = threshs(4);
cfg_CARACAS.thresh_Rampl = threshs(5);
cfg_CARACAS.thresh_bpm = threshs(6:7);

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
Weighted_accuracy = 0.25*Sensitivity + 0.75*Specificity;
% Accuracy = ( sum(toplot == 1) + sum(toplot == 3) ) / (sum(toplot == 1) + sum(toplot == 2) + sum(toplot == 3) + sum(toplot == 4) );

myerr = 1 - Weighted_accuracy;
end

