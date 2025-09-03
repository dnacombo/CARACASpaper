setpath_ds003690

compare_truth = 'MANUAL';
compare_with = 'CARACAS'; %'ICLabel'; % 'CARACAS'; % 'CORR';

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


%% Number of folds
K = 5;
cfg_SASICA = SASICA('getdefs');
cfg_SASICA.CARACAS.enable = 1;
cfg_CARACAS = cfg_SASICA.CARACAS;
param_vec = [cfg_CARACAS.thresh_sk cfg_CARACAS.thresh_ku cfg_CARACAS.thresh_PQ cfg_CARACAS.thresh_RR cfg_CARACAS.thresh_Rampl cfg_CARACAS.thresh_bpm];

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
param_names = {'thresh_sk', 'thresh_ku1', 'thresh_PQ', 'thresh_RR', 'thresh_Rampl', 'thresh_bpm1', 'thresh_bpm2', 'TestError'};
results_table = [bestparams, testerr'];
mean_row = [mean(bestparams), mean(testerr)];

% Create table
T = array2table([results_table; mean_row], 'VariableNames', param_names);
T.Properties.RowNames = [arrayfun(@(x) sprintf('Fold_%d', x), 1:K, 'UniformOutput', false), {'Mean'}];

% Display table
disp('Cross-validation Results:')
disp(T)

%%

function myerr = CARACASerr(threshs, fs, compare_truth, compare_with)

cfg_SASICA = SASICA('getdefs');
cfg_CARACAS = cfg_SASICA.CARACAS;
cfg_CARACAS.thresh_sk = threshs(1);
cfg_CARACAS.thresh_ku = threshs(2);
cfg_CARACAS.thresh_PQ = threshs(3);
cfg_CARACAS.thresh_RR = threshs(4);
cfg_CARACAS.thresh_Rampl = threshs(5);
cfg_CARACAS.thresh_bpm = threshs(6:7);

CORRrej = zeros(numel(fs),numel(fs(1).CORR.rej));CARACASrej = zeros(numel(fs),numel(fs(1).CARACAS.rej));SASICARACASrej = zeros(numel(fs),numel(fs(1).SASICARACAS.rej));
for i = 1:numel(fs)
    CORRrej(i,:) = fs(i).CORR.rej;

    CARACASrej(i,:) = CARACAS_rethresh(CARACASrej(i,:), fs(i).CARACAS, cfg_CARACAS);

    SASICARACASrej(i,:) = CARACAS_rethresh(SASICARACASrej(i,:), fs(i).SASICARACAS, cfg_CARACAS);
    MANUALrej(i,:) = fs(i).MANUAL.rej;
end
switch compare_truth
    case 'MANUAL'
        truthrej = MANUALrej;
    case 'SASICARACAS'
        truthrej = SASICARACASrej;
    case 'CARACAS'
        truthrej = CARACASrej;
    case 'CORR'
        truthrej = CORRrej;
end
switch compare_with
    case 'MANUAL'
        withrej = MANUALrej;
    case 'SASICARACAS'
        withrej = SASICARACASrej;
    case 'CARACAS'
        withrej = CARACASrej;
    case 'CORR'
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

