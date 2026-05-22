dirroot = fileparts(fileparts(mfilename('fullpath')));
dircode = fullfile(dirroot,'code');
dirbids = fullfile(dirroot,'ds003690');
if not(exist(dirbids,'file'))
    error(['Download dataset from https://github.com/OpenNeuroDatasets/ds003690.git and put it in ' dirbids])
end
dirderiv = fullfile(dirbids,'derivatives');


dirout = fullfile(dirderiv,['CardiClassif' which_data]);

addpath(fullfile(dircode,'MiscMatlab'))
addpath(fullfile(dircode,'MiscMatlab', 'stats'))
addpath(fullfile(dircode,'MiscMatlab','plot'))
addpath(fullfile(dircode,'MiscMatlab','plot','panel'))

addpath(fullfile(dircode,'SASICA'))
rm_frompath('eeglab') % avoid interference with other versions
rm_frompath('fieldtrip')
% add local fieldtrip install
addpath(fullfile(dircode,'fieldtrip'))
addpath(fullfile(dircode,'fieldtrip', 'external','eeglab'))
ft_defaults

addpath(fullfile(dircode,'SASICA','CARACAS'));
addpath(fullfile(dircode,'SASICA','CARACAS', 'heart_functions'));

cfg_SASICA = SASICA('getdefs');
addpath(fullfile(dircode, "SASICA",'eeglab'))
addpath(genpath(fullfile(dircode, "SASICA",'eeglab', 'functions')))
addpath(fullfile(dircode, "SASICA",'eeglab', 'plugins', 'ICLabel'))
addpath(fullfile(dircode, "SASICA",'eeglab', 'plugins', 'firfilt'))

mymkdir(dirout)
