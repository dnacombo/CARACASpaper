%
%PAL_hdi  Find high-density interval (aka credible interval)
%
%   syntax: [hdi] = PAL_hdi(x, pdf, width)
%
%   Find high-density interval (HDI, aka credible interval) in probability 
%   mass function or discretized probablity density function pdf defined 
%   over quantitative variable x.
%
%Input:
%
%   'x': Vector containing discrete values of the random variable.
%
%   'pdf': probability mass function over variable x
%
%   'width': width of high density interval reported. Follow by scalar. If
%       supplied value is in (0, 1) it is interpreted as a proportion, if
%       it is in [1, 100) it will be interpreted as percentage. Default:
%       0.95.
%
%Output:
%
%   'hdi': high-density interval as n x 2 array. Typically, n will equal 1
%       and the 1 x 2 vector will contain lower and upper limit of the 
%       single interval that corresponds to HDI. In case HDI consists of 
%       multiple regions, the array will have as many rows as there are
%       separate regions (see example 2).
%
%Example 1:
%
%   x = randn(1,10000); %draw sample from standard normal distribution
%   [grid,pdf] = PAL_kde(x);  %estimate density function
%   hdi = PAL_hdi(grid,pdf,68.28)
%   plot(grid,pdf);
%   line([hdi(1) hdi(2)],[pdf(find(grid == hdi(1))), pdf(find(grid == hdi(1)))]);
%
%Example 2:
%
%   x = [randn(1,10000),randn(1,7000)+5]; %draw samples from two normal distributions
%   [grid,pdf] = PAL_kde(x);  %estimate density function
%   hdi = PAL_hdi(grid,pdf,.7)
%   plot(grid,pdf);
%   height = pdf(find(grid == hdi(2)));
%   line([hdi(1,1) hdi(1,2)],[height, height]);
%   line([hdi(2,1) hdi(2,2)],[height, height]);
%
%Introduced: Palamedes version 1.10.10 (NP)

function hdi = PAL_hdi(x,pdf,width)

if width >= 1
    width = width/100;
end

[dis I] = sort(pdf/sum(pdf),'descend');
dis = cumsum(dis);
dis = sort(I(dis < width));
discont = find((dis(2:end) - dis(1:end-1))~=1);

begins = [1 discont+1];
ends = [discont length(dis)];

for region = 1:length(begins)
    hdi(region,:) = [x(dis(begins(region))),x(dis(ends(region)))];
end