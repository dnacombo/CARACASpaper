dirroot = '/network/iss/cenir/analyse/meeg/CARACAS/Test_Max/';
dircode = fullfile(dirroot,'code');
dirbids       = fullfile(dirroot,'ds003690');
dirderiv = fullfile(dirbids,'derivatives');
dirout = fullfile(dirderiv,'CardiClassif');


addpath(fullfile(dircode,'MiscMatlab'))
addpath(fullfile(dircode,'MiscMatlab', 'stats'))
addpath(fullfile(dircode,'MiscMatlab/plot/'))
addpath(fullfile(dircode,'MiscMatlab/plot/panel'))

addpath(fullfile(dircode,'Palamedes'))
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
addpath(fullfile(dircode,'Cardiac_IC_labelling', 'heart_functions/'));

cfg_SASICA = SASICA('getdefs');
addpath(genpath(fullfile(dircode, "SASICA",'eeglab')))

