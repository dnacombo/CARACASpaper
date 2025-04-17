%
%PAL_pdfExGaussian  Ex-Gaussian probability density
%
%syntax: y = PAL_pdfExGaussian(x, mu, sigma, tau)
%
%Returns the probability density of the ExGaussian distribution with 
%   parameters 'mu', 'sigma', and 'tau' evaluated at 'x'. 
% 
%   'x' may be array of any size 
% 
%   'mu', 'sigma', and 'tau' should either be scalars or arrays in size 
%       equal to 'x'.
%
%   'mu' and 'sigma' correspond to the mean and standard deviation of the
%       Gaussian component. 'tau' corresponds to the 'relaxation time' of
%       the exponential component.
%
%Examples:
%   
%   x = 0:200;
%
%   y = PAL_pdfExGaussian(x,30,10,30); 
%
%Introduced: Palamedes version 1.11.11 (NP)

function y = PAL_pdfExGaussian(x,mu,sigma,tau)

y = exp((-1.*(x-mu)./tau) + (sigma./tau).^2).*PAL_CumulativeNormal([0 1],(x-mu)./sigma - (sigma/tau));
