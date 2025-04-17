%
%PAL_kde2  Bivariate Kernel Density Estimation
%
%   syntax: [gridx, gridy, pdf, cdf] = PAL_kde2(x,y,{optional boundaries})
%
%   Estimates joint probability density for bivariate distribution based on
%   randomly sampled, paired values in vectors x and y.
%
%Input:
%
%   'x': Vector containing sampled values of one of the two random 
%       variables whose joint density is to be estimated.
%
%   'y': Vector containing sampled values of the other random variable 
%       whose joint density is to be estimated. Values in x and y must 
%       represent pairs of sampled values.
%
%Output:
%
%   'gridx': Values of random variable x
%
%   'gridy': Values of random variable y
%
%   'pdf': Estimated probability density function across 'gridx' and 
%       'gridy'
%
%   'cdf': Estimated cumulative distribution function across 'gridx' and 
%       'gridy'
%   
%   If underlying random variables are constrained to a finite interval 
%   (e.g., probabilities to [0 1]) user may supply a 2x2 array containing 
%   lower and upper boundaries for 'x' in the first row and those for 'y' 
%   in the second row. Use of e.g., [-Inf, 10; 0 Inf] is allowed. Using 
%   empty array will result in the default [-Inf Inf; -Inf Inf].
%
%   Calculations are somewhat intensive and can noticeably slow down
%   execution of a routine if this function is called. For a quicker, less 
%   accurate probability density supply optional argument: 'quick'
%
%Example 1:
%
%   xy = PAL_rndBiNormal(10000,[1 2],[4 6; 6 16]); %draw samples from 
%       %bivariate normal distribution
%   [gridx,gridy,pdf,cdf] = PAL_kde2(xy(:,1),xy(:,2));
%   [gridx,gridy] = ndgrid(gridx,gridy);
%   surf(gridx,gridy,pdf);
%
%Example 2:
%
%   x = rand(1,100000);  %draw samples from bivariate uniform distribution
%   y = rand(1,100000);
%   [gridx,gridy,pdf] = PAL_kde2(x,y,[0 1;0 1]);
%   [gridx,gridy] = ndgrid(gridx,gridy);
%   surf(gridx,gridy,pdf);
%   set(gca,'zlim',[0 max(max(pdf))]);
%
%Example 3:
%
%   x = rand(1,100000);  %draw samples from bivariate uniform distribution
%   y = rand(1,100000);
%   [gridx,gridy,pdf] = PAL_kde2(x,y,[-Inf 1;0 Inf]); %weird example but 
%                                                  %demonstrates boundaries
%   [gridx,gridy] = ndgrid(gridx,gridy);
%   surf(gridx,gridy,pdf);
%
%Introduced: Palamedes version 1.11.7 (NP)
%Modified: Palamedes version 1.11.7 (See History.m)


function [centersx, centersy, pdf, cdf] = PAL_kde2(x, y,varargin)

bound = [-Inf Inf; -Inf Inf];
nbins = 101;

if ~isempty(varargin)
    bound = varargin{1};
    if isempty(bound)
        bound = [-Inf Inf;-Inf Inf];
    end
    if length(varargin) == 2
        if strncmpi(varargin{2}, 'quick',1)
            nbins = 31;
        end
    end
end

binwidth = [(max(x)-min(x))/(nbins-1),(max(y)-min(y))/(nbins-1)];
h = [std(x)*length(x).^(-1/6), std(y)*length(y).^(-1/6)]; %Kernel bandwidth suggested by Silverman's (1986) rule of thumb 
                                                          %(https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation, retrieved 2/5/2023))
kernelpad = 4*h;
padding = ceil(kernelpad./binwidth);

edgesx = min(x)-padding(1)*binwidth(1):binwidth(1):max(x)+padding(1)*binwidth(1);
edgesy = min(y)-padding(2)*binwidth(2):binwidth(2):max(y)+padding(2)*binwidth(2);
centersx = edgesx(1:end-1)+binwidth(1)/2;
centersy = edgesy(1:end-1)+binwidth(2)/2;
[gridx, gridy] = ndgrid(centersx,centersy);

n = histcounts2(x,y,edgesx,edgesy);

pdf = zeros(size(gridx));

for xI = 1:size(gridx,1)
    for yI = 1:size(gridx,2)
        pdf = pdf+PAL_pdfBiNormal(gridx,gridy,[gridx(xI,yI),gridy(xI,yI)],[h(1), h(2)]).*n(xI,yI);
    end
end

if ~isinf(bound(1,1))
    I = find(centersx>=bound(1,1),1,'first');
    pdf(I:2*I-2,:) = pdf(I:2*I-2,:)+pdf(I-1:-1:1,:);
    pdf = pdf(I:end,:);
    centersx = centersx(I:end);
end
if ~isinf(bound(1,2))    
    I = find(centersx<=bound(1,2),1,'last');
    pdf(2*I-end+1:I,:) = pdf(2*I-end+1:I,:)+pdf(end:-1:I+1,:);
    pdf = pdf(1:I,:);
    centersx = centersx(1:I);
end
if ~isinf(bound(2,1))
    I = find(centersy>=bound(2,1),1,'first');
    pdf(:,I:2*I-2) = pdf(:,I:2*I-2)+pdf(:,I-1:-1:1);
    pdf = pdf(:,I:end);
    centersy = centersy(I:end);
end
if ~isinf(bound(2,2))
    I = find(centersy<=bound(2,2),1,'last');
    pdf(:,2*I-end+1:I) = pdf(:,2*I-end+1:I)+pdf(:,end:-1:I+1);
    pdf = pdf(:,1:I);
    centersy = centersy(1:I);
end

pdf = pdf/sum(sum(n));
cdf = cumsum(cumsum(pdf),2).*binwidth(1).*binwidth(2);