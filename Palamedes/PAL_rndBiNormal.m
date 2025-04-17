%
%PAL_rndBiNormal  Random draws from Bivariate Normal distribution
%
%syntax: x = PAL_rndBiNormal(n, Mean, COV)
%        or:
%        x = PAL_rndBiNormal(n, Mean, SD, r)
%
%Returns n random draws from bivariate normal distribution with specified
%   parameters.
%
%Input:
%
%   'n': positive integer specifying the number of draws to be taken.
% 
%   'Mean' (optional): 1x2 vector specifying the means of the two variables
%
%   'COV' (optional): 2x2 array specifying the covariance matrix
%
%   'SD' (optional): 1x2 vector specifying the standard deviations of the
%       two variables
%
%   'r' (optional): scalar specifying the correlation between the two 
%       variables
%
%   Only 'n' must be supplied, a single additional argument will be
%       interpreted as 'Mean' above, a second additional argument will be
%       interpreted as 'COV' or 'SD' above (depending on size), a third
%       additional argument will be interpreted as 'r' above. If not
%       supplied, 'Mean' will be assumed to be [0 0], 'COV' will be assumed
%       to be 2x2 identity matrix (equivalently, 'SD' will be assumed to be
%       [1 1] and 'r' will be assumed to be 0).
%
%Output:
%
%   'x': n x 2 array of paired values drawn from bivariate normal
%       distribution.
%
%Example 1:
%
%   The following three calls are equivalent:
%
%   x = PAL_rndBiNormal(10)
%
%   x = PAL_rndBiNormal(10,[0 0],[1 1],0)
%
%   x = PAL_rndBiNormal(10,[0 0],[1 0; 0 1])
%
%Wxample 2:
%
%   Given, let's say:
%
%   Mean = [1 2];
%   SD = [3 4];
%   r = 0.5;
% 
%   COV = [SD(1)^2, r*SD(1)*SD(2); r*SD(1)*SD(2) SD(2)^2];
% 
%   the following two calls are equivalent:
% 
%   x = PAL_rndBiNormal(100, Mean, SD, r)
%
%   x = PAL_rndBiNormal(100, Mean, COV)
%
%   scatter(x(:,1), x(:,2)) %display draws
%
%Introduced: Palamedes version 1.11.7 (NP)

function [x] = PAL_rndBiNormal(n,varargin)

mu = [0 0];
sd = [1 1];
r = 0;

if ~isempty(varargin)
    NumOpts = length(varargin);
    mu = varargin{1};
    if NumOpts > 1        
        if isvector(varargin{2})
            sd = varargin{2};
            if NumOpts > 2
                r = varargin{3};
            end
        else
            cov = varargin{2};
            sd = [sqrt(cov(1,1)), sqrt(cov(2,2))];
            r = cov(2,1)/(sd(1)*sd(2));
        end
    end
end

x = randn(n,1)*sd(1) + mu(1);
x(:,2) = randn(n,1).*sqrt(((1 - r.^2).*sd(2).^2))+ mu(2) + (sd(2)/sd(1)).*r.*(x(:,1) - mu(1)); 
