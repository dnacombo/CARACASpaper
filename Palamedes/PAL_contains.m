%
%PAL_contains  Emulates (some) functionality of Matlab's 'contains'
%   function
%
%   For compatibility with Octave
%
%   syntax: [patternFound] = PAL_contains(str, pattern)
%
%   PAL_contains returns logical true if string 'pattern' is contained in 
%       string 'str' or logical false otherwise.
%
% Introduced: Palamedes version 1.10.0 (NP)
% Modified: Palamedes version 1.11.8 (See History.m)

function [patternFound] = PAL_contains(str,pattern)

if iscellstr(str)
    for I = 1:length(str)        
        patternFound(I) = ~isempty(strfind(str{I},pattern));
    end
else
    patternFound = ~isempty(strfind(str,pattern));
end