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
EEG.icaweights = comp.unmixing;
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
EEG.chaninfo.nosedir = '+X';

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
EEG.icachansind = 1:size(EEG.icaact,1);

if isfield(cfg,'reject')
    EEG.reject = cfg.reject;
end
end