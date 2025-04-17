%
%PAL_pdfBiNormal  Bivariate Normal probability density
%
%syntax: z = PAL_pdfBiNormal(x, y, Mean, COV)
%        or:
%        z = PAL_pdfBiNormal(x, y, Mean, SD, r)
%
%Returns the probability density of the bivariate normal distribution with 
%   specified parameters evaluated at paired values x and y.
%
%Input:
%
%   'x': array of any size specifying the values of one of the two 
%       variables at which the density should be evaluated.
% 
%   'y': array of same size as x specifying the values of the second 
%       variable at which the density should be evaluated.
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
%   Only 'x' and 'y' must be supplied, a single additional argument will 
%       be interpreted as 'Mean' (above), a second additional argument will
%       be interpreted as 'COV' or 'SD' above (depending on size), a third
%       additional argument will be interpreted as 'r' above. If not
%       supplied, 'Mean' will be assumed to be [0 0], 'COV' will be assumed
%       to be 2x2 identity matrix (equivalently, 'SD' will be assumed to be
%       [1 1] and 'r' will be assumed to be 0).
%
%Output:
%
%   'z': probability density at x, y
%
%Example 1:
%
%   [xgrid ygrid] = ndgrid(-10:.1:10,-10:.1:10);
%   z = PAL_pdfBiNormal(xgrid,ygrid,[2 -2],[2 1],-.3);
%   surf(xgrid,ygrid,z);
%
%   The following three calls are equivalent:
%
%   z = PAL_pdfBiNormal([1 2])
%
%   z = PAL_pdfBiNormal([1 2],[0 0],[1 1],0)
%
%   z = PAL_pdfBiNormal([1 2],[0 0],[1 0; 0 1])
%
%Example 2:
%
%   Given, let's say:
%
%   X = [2 3];
%   Mean = [1 2];
%   SD = [3 4];
%   r = 0.5;
%
%   COV = [SD(1)^2, r*SD(1)*SD(2); r*SD(1)*SD(2) SD(2)^2];
%
%   the following two calls are equivalent:
%
%   z = PAL_pdfBiNormal(X(1), X(2), Mean, SD, r)
%
%   z = PAL_pdfBiNormal(X(1), X(2), Mean, COV)
%
%Introduced: Palamedes version 1.11.7 (NP)

function z = PAL_pdfBiNormal(x, y, varargin)

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

x = (x-mu(1))./sd(1);
y = (y-mu(2))./sd(2);
z = 1./(2*pi.*sd(1).*sd(2).*sqrt(1-r.^2)).*exp((-1./(2.*(1-r.^2))).*(x.^2-2.*r.*x.*y+y.^2));