%
%PAL_pdfBeta  Beta probability density
%
%syntax: y = PAL_pdfBeta(x, a, b, {optional argument})
%
%Returns the probability density of the beta distribution with parameters 
%   'a' and 'b' evaluated at 'x'.
%
%   'x' may be array of any size 
% 
%   'a' and 'b' should either be scalars or arrays in size equal to 'x'.
%
%   By default, 'a' and 'b' are interpreted as the 'alpha' and 'beta'
%       parameters in the standard formulation of the beta distribution.
%       Same can be accomplished by providing optional argument 'ab'.
%       If the optional argument 'meanandconcentration' (or just 'mean') is 
%       provided, 'a' and 'b' are interpreted as the mean and concentration 
%       parameters of the beta distribution. If the optional argument
%       'modeandconcentration' (or just 'mode') is provided, 'a' and 'b' 
%       are interpreted as the mode and concentration parameters of the 
%       beta distribution. The latter two cases are perhaps easier to 
%       intuit, as concentration can be interpreted as roughly equal to 
%       'sample size'. For example, a beta distribution with mean parameter 
%       equal to 0.03 and concentration parameter equal to 100 is 
%       approximately as informative as having observed 3 'successes' out 
%       of 100 bernouilli trials.
%
%Note that, for some parameter values, density may evaluate to Inf at x = 0
%   and/or x = 1. It may be desirable to avoid this.
%
%Examples:
%   
%   x_grain = 0.01;
%   x = x_grain/2:x_grain:1-x_grain/2; %Avoid x = 0 and x = 1
%
%   y = PAL_pdfBeta(x,3,97); %Gives beta with alpha = 3 and beta = 97
%
%   y = PAL_pdfBeta(x,0.03,100,'mean'); %Gives beta with 
%       %mean = 0.03 and concentration = 100 (which is equivalent to
%       %previous example)
%
%   y = PAL_pdfBeta(x,0.3,50,'mode'); %Gives beta with 
%       %mode = 0.2 and concentration = 50
%
%Introduced: Palamedes version 1.11.11 (NP)

function y = PAL_pdfBeta(x,p1,p2,varargin)

a = p1;
b = p2;

valid = 0;
if ~isempty(varargin)
    if strncmpi(varargin{1},'meanandconcentration',2)
        a = p1.*p2;
        b = (1-p1).*p2;
        valid = 1;
    end
    if strncmpi(varargin{1},'modeandconcentration',2)
        a = p1.*(p2-2)+1;
        b = (1-p1).*(p2-2)+1;
        valid = 1;
    end
    if strncmpi(varargin{1},'ab',2)
        valid = 1;
    end
    if ~valid
        warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{1});
    end
end

y = (gamma(a+b)./(gamma(a).*gamma(b))).*x.^(a-1).*(1-x).^(b-1);