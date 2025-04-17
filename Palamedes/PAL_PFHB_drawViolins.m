%
%PAL_PFHB_drawViolins  Display violin plots for (some) parameters in 
%   analysis performed by PAL_PFHB_fitModel
%
%   syntax: PAL_PFHB_drawViolins(pfhb, {optional: parameter},{more optional arguments})
%
%   Dot corresponds to central tendency ({mode}, median, or mean) in 
%   marginal posterior distribution, line covers 68% high-density interval 
%   (hdi, or credible interval), curves show posterior across its 95% hdi.
%
%Input: 
%
%   'pfhb': Analysis structure created by PAL_PFHB_fitModel
%
%   PAL_PFHB_drawViolins also accepts optional arguments    
%
%   If no parameter is specified (se below), this routine will display 
%       violin plots for what are likely the parameters of most interest to 
%       you:
%
%       If single-subject data, violins for location ('threshold'), slope,
%       guess and lapse parameters are shown (if parameter was free) for
%       all conditions/effects in analysis.
%   
%       If multiple-subject data, violins for the hyperparameters of
%       location ('threshold'), slope, guess and lapse parameters are shown 
%       (if parameter was free) for all conditions/effects in analysis.
%
%Optional arguments:
%
%   centralTendency: followed by additional argument 'mean', {'mode'},
%       'median', uses indicated measure of central tendency for the dot in 
%       the plot.
%
%   To inspect any parameter, provide the parameter as the first optional 
%       argument (i.e., immediately following the mandatory argument 
%       specifying the analysis structure). To see a listing of free 
%       parameters in the model use PAL_PFHB_drawViolins(pfhb,'list');
%
%       If the first optional argument is used to specify a parameter,
%       results can be limited to a subset of conditions/effects and/or
%       subjects. To specify a subset of conditions or effects, use a pair
%       of arguments, the first of which is either 'conditions' or
%       'effects' (for this purpose, these terms are used interchangeably),
%       followed by a vector specifying the subset (e.g.,
%
%       PAL_PFHB_drawViolins(pfhb,'a','conditions', [1 3:5]);
%
%       In order to limit results to a subset of subjects use argument
%       'subjects', followed by a vector specifying the subset.
%       Arguments 'conditions' (or 'effects') and 'subjects' can be used
%       simultaneously.
%
%For examples of use see any of the PAL_PFHB demos in the PalamedesDemos
%   folder or visit:
%   www.apalamedestoolbox.org/hierarchicalbayesian.html
%
%Introduced: Palamedes version 1.10.10 (NP)
%Modified: Palamedes version 1.11.8, 1.11.9 (see History.m)

function [] = PAL_PFHB_drawViolins(pfhb,varargin)

paramSupplied = false;
paramsModeled = unique(pfhb.summStats.linList.p(~strcmp(pfhb.summStats.linList.p,'deviance')));
centTend = 'mode';

if ~isempty(varargin)
    if strcmp(varargin{1},'list')
        disp(pfhb.model.paramsList);
        return;
    else        
        if ~any(strcmp(varargin{1},paramsModeled))
            if strncmpi(varargin{1},'central tendency',4)
                centTend = varargin{2};
            else
                message = ['Parameter ',param,' was not modeled in this analysis. Use one of the parameters listed below (e.g., ''',char(pfhb.summStats.linList.p{1}) ,'''). ',pfhb.model.paramsList];
                error('PALAMEDES:invalidOption',message);   
            end
            if any(strncmpi(varargin,'subjects',4)) || any(strncmpi(varargin,'conditions',4)) || any(strncmpi(varargin,'effects',4))
                warning('PALAMEDES:invalidOption','''subjects'' and ''conditions/effects'' are not valid options unless a parameter is specified. Ignored.');
            end
        else
            param = varargin{1};
            paramSupplied = true;
        end
    end
end

if ~isfield(pfhb.summStats.(pfhb.summStats.linList.p{1}),'hdi68')    %pfhb created by Palamedes < version 1.10.10?    
    pfhb.summStats = PAL_PFHB_updateOldStyleStructure(pfhb.summStats);
end

if ~paramSupplied

    params = pfhb.model.parameters(find(~PAL_contains(pfhb.model.parameters,'_actual')));
    if pfhb.model.Nsubj > 1
        params = params(find(PAL_contains(params,'mu')));
    end
    if pfhb.model.Ncond > 1
        params = fliplr(params);
    end

    for paramIndex = 1:length(params);
        param = params{paramIndex};
        numEffects(paramIndex) = pfhb.model.(param(1)).Nc;
        if PAL_mmType(pfhb.model.(param(1)).c) == 1;
            xLabel{paramIndex} = 'condition';
        else
            xLabel{paramIndex} = 'effect';
        end
        switch param
            case 'a'
                yLabel{paramIndex} = 'location (''threshold'')';
            case 'b'
                yLabel{paramIndex} = 'slope';
            case 'g'
                yLabel{paramIndex} = 'guess';
            case 'l'
                yLabel{paramIndex} = 'lapse';
            case 'amu'
                yLabel{paramIndex} = 'location mean';
            case 'bmu'
                yLabel{paramIndex} = 'slope mean';
            case 'gmu'
                yLabel{paramIndex} = 'guess mean';
            case 'lmu'
                yLabel{paramIndex} = 'lapse mean';
        end
    end
    maxEffects = max(numEffects);

    if pfhb.model.Ncond == 1
        figpos = [100 100 length(params)*200+50 500+20];
        axesadv = [200 0];
        axessize = [125 425];
    else           
        figpos = [100 100 maxEffects*60+100 length(params)*200+20];
        axesadv = [0 190];
        axessize = [maxEffects*60 140];
    end

    figure('units','pixels','position',figpos,'color','w');

    subject = 1;

    S = warning('off','PALAMEDES:invalidOption');
    for paramIndex = 1:length(params)
        param = params{paramIndex};

            axes('units','pixels','position',[[75 50]+((paramIndex-1)*axesadv) axessize]);
            set(gca,'fontsize',12);
            set(gca,'xlim',[.5 maxEffects+.5]);
            set(gca,'xtick',[1:numEffects(paramIndex)]);
            xlabel(xLabel{paramIndex});
            ylabel(yLabel{paramIndex});

            hold on            

            for effect = 1:size(pfhb.summStats.(param).mode,1)
                plot(effect,pfhb.summStats.(param).(centTend)(effect,subject),'o','color','k','markerfacecolor','k','markersize',10);
                line([effect effect],[pfhb.summStats.(param).hdi68(effect,subject,1) pfhb.summStats.(param).hdi68(effect,subject,2)],'color','k','linewidth',2);
                [stats samples] = PAL_PFHB_inspectParam(pfhb,param,'effect',effect,'subject',subject,'nofig','hdi',95);
                s = samples(:);  
                if param(1) == 'g' || param(1) == 'l'
                    bounds = [0 1];
                else
                    bounds = [-Inf Inf];
                end
                [grid pdf] = PAL_kde(s,bounds);                                    
                pdf = pdf(grid > stats.hdi(1,1) & grid < stats.hdi(end,2));
                grid = grid(grid > stats.hdi(1,1) & grid < stats.hdi(end,2));
                pdf = .4*pdf/max(pdf);
                plot(effect-pdf,grid,'color','k','linewidth',2);
                plot(effect+pdf,grid,'color','k','linewidth',2);   
            end
            if paramIndex == 1
                ylim = get(gca,'ylim');
            end
    end
    
else
    
    fullHyperParamList = {'amu','amu_actual','asigma','bmu','bmu_actual','bsigma','gmu','gmu_actual','gkappa','lmu','lmu_actual','lkappa'};

    if any(strcmp(param,fullHyperParamList))
        hyper = true;
        subjectsList = 1;
    else
        hyper = false;
        subjectsList = 1:pfhb.model.Nsubj;
    end

    if ~PAL_contains(param,'_actual')
        effectsList = 1:pfhb.model.(param(1)).Nc;    
    else
        effectsList = 1:pfhb.model.Ncond; 
    end

    subjectsSupplied = false;
    effectsSupplied = false;
    
    if length(varargin) > 1
        NumOpts = length(varargin);
        n = 2;
        while n <= NumOpts
            valid = 0;
            if strncmpi(varargin{n},'central tendency',4)
                if any(strcmpi(varargin{n+1},{'mode','mean','median'}))
                    centTend = varargin{n+1};
                else
                    warning('PALAMEDES:invalidOption','%s is not a valid option for a central tendency. Ignored.',varargin{n+1});
                end
                valid = 1;
            end
            if strncmpi(varargin{n},'subjects',4)
                if ~hyper
                    subjectsList = varargin{n+1};            
                    subjectsSupplied = true;                    
                else
                    warning('PALAMEDES:invalidOption','''subjects'' is not a valid option for a hyper parameter. Ignored.');
                end
                valid = 1;
            end
            if any(strncmpi(varargin{n},{'conditions','effects'},4))
                effectsList = varargin{n+1};      
                effectsSupplied = true;
                valid = 1;
            end                        
            if valid == 0
                warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
            end
            n = n + 2;
        end
    end

    S = warning('off','PALAMEDES:invalidOption');

    if length(subjectsList) > 6 && length(effectsList) > 1
        message = ['When displaying results for multiple conditions or effects, I can only display 6 subjects at most.',char(10), ...
            'Use ''subjects'' option to display results for other subjects or use ''conditions'' option to display results for a single condition.'];
        disp(message);
        subjectsList = subjectsList(1:6);
    end

    Nsubjects = length(subjectsList);
    Neffects = length(effectsList);

    if Neffects == 1
        axessize = [Nsubjects*60 120]
        xlim = [0 Nsubjects+1];
        xtick = 1:Nsubjects;
        if ~hyper                                               %Fixed here
            xticklabel = strsplit(int2str(subjectsList));
        else
            xticklabel = strsplit(int2str(effectsList));
        end
        if hyper
            if PAL_mmType(pfhb.model.(param(1)).c) == 1 || PAL_contains(param,'_actual')
                xlabeltext = 'condition';
            else
                xlabeltext = 'effect';
            end
        else
            xlabeltext = 'subject';
        end
        effects = effectsList.*ones(1,length(subjectsList));    %Fixed here used to be: effects = ones(1,length(subjectsList));
        subjects = subjectsList;
        Naxes = 1;
    else
        axessize = [Neffects*60 120];
        xlim = [0 Neffects+1];
        xtick = 1:Neffects;
        xticklabel = strsplit(int2str(effectsList));
        if PAL_mmType(pfhb.model.(param(1)).c) == 1 || PAL_contains(param,'_actual') 
            xlabeltext = 'condition';                                                   
        else
            xlabeltext = 'effect';
        end                                                                     
        Naxes = Nsubjects;
        if Naxes > 1
            subjectsList = fliplr(subjectsList);
        end
        effects = repmat(effectsList,[length(subjectsList) 1]);
        subjects = repmat(subjectsList',[1 length(effectsList)]);

    end

    figure('units','pixels','position',[100 100 axessize(1)+140 Naxes*150+100],'color','w');

    for axIndex = 1:Naxes
        ax(axIndex) = axes('units','pixels','position',[[75 50]+((axIndex-1)*[0 150]) axessize]);    
        set(ax(axIndex),'xlim',xlim,'xtick',xtick,'xticklabel',xticklabel);
        xlabel(xlabeltext);
        if Naxes > 1
            ylabel(['subject ',int2str(subjectsList(axIndex))]);
        end
        hold on;
        for xIndex = 1:size(effects,2)
            plot(ax(axIndex),xIndex,pfhb.summStats.(param).(centTend)(effects(axIndex,xIndex),subjects(axIndex,xIndex)),'o','color','k','markerfacecolor','k','markersize',6);
            line(ax(axIndex),[xIndex xIndex],[pfhb.summStats.(param).hdi68(effects(axIndex,xIndex),subjects(axIndex,xIndex),1) pfhb.summStats.(param).hdi68(effects(axIndex,xIndex),subjects(axIndex,xIndex),2)],'color','k','linewidth',2);
            [stats samples] = PAL_PFHB_inspectParam(pfhb,param,'effect',effects(axIndex,xIndex),'subject',subjects(axIndex,xIndex),'nofig','hdi',95);
            s = samples(:);  
            bounds = [-Inf Inf];
            if any(strcmp(param,{'g','l','gmu','lmu'}))
                bounds = [0 1];
            end
            if any(strcmp(param,{'asigma','bsigma','gkappa','lkappa'}))
                bounds = [0 Inf];
            end
            [grid pdf] = PAL_kde(s,bounds);        
            pdf = pdf(grid > stats.hdi(1,1) & grid < stats.hdi(end,2));
            grid = grid(grid > stats.hdi(1,1) & grid < stats.hdi(end,2));
            pdf = .4*pdf/max(pdf);
            plot(ax(axIndex),xIndex-pdf,grid,'color','k','linewidth',2);
            plot(ax(axIndex),xIndex+pdf,grid,'color','k','linewidth',2);   
        end
        if axIndex == Naxes
            ylim = get(ax(axIndex),'ylim');
            text(mean([xlim(1),xlim(2)]),ylim(2)+(ylim(2)-ylim(1))/10,['parameter: ',param],'horizontalalignment','center','verticalalignment','bottom','fontsize',12,'interpreter','none')
        end
    end
end    
S = warning('on','PALAMEDES:invalidOption');