%
%PAL_PFHB_updateOldStyleStructure  Update old style (before version 
%   1.10.10) PFHB analysis structure to be compatible with
%   PAL_PFHB_drawViolins.
%
%   Analysis structures (e.g., 'pfhb') created by older versions of 
%   PAL_PFHB_fitModel had entries pfhb.summStats.[parameter].HDI68low, 
%   ...HDI68high, ...HDI95low, and ...HDI95high. Currently, analysis
%   structures store HDIs in entries pfhb.summStats.[parameter].hdi68 and
%   ...hdi95. This function simply reads old style HDIs and creates the new
%   style entries.
%
%Introduced: Palamedes version 1.10.10 (NP)

function [summStats] = PAL_PFHB_updateOldStyleStructure(summStats)

if ~isfield(summStats.(summStats.linList.p{1}),'hdi68')
        for index = 1:length(summStats.linList.p)        
        summStats.(summStats.linList.p{index}).hdi68(summStats.linList.c(index),summStats.linList.s(index),:) = [summStats.(summStats.linList.p{index}).HDI68low(summStats.linList.c(index),summStats.linList.s(index)) summStats.(summStats.linList.p{index}).HDI68high(summStats.linList.c(index),summStats.linList.s(index))];
        summStats.(summStats.linList.p{index}).hdi95(summStats.linList.c(index),summStats.linList.s(index),:) = [summStats.(summStats.linList.p{index}).HDI95low(summStats.linList.c(index),summStats.linList.s(index)) summStats.(summStats.linList.p{index}).HDI95high(summStats.linList.c(index),summStats.linList.s(index))];
    end
end