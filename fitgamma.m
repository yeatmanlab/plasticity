function coefs = fitgamma(x,y,init,show)
% Fits a Poisson curve to data points
%
% coefs = fitPoissonCurve(x,y,show)
%
% Nonlinear search, with multiple starting points is used to fit a Poisson
% curve to lifespan data. Equation for a Poisson curve is c+a*x*e^-bx 
%
% Inputs:
% x and y are vectors of data. Show is a logical denoting whether to plot
% the fit.
%
% Outputs:
% coefs(1) = a
% coefs(2) = b
% coefs(3) = c
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.

if ~exist('init')
    init = [1 2 1];;
end
pfun = @(p,x) evalgamma(p,x)
options = optimset('Display','off');
for ii = 1:size(init,1)
    [coefs(ii,:),resnorm(ii,:)] = lsqcurvefit(pfun,init,x,y,[],[],options);
end
% Find best solution
[~,ind]=min(resnorm);
coefs=coefs(ind,:);
if exist('show','var') && ~isempty(show) && show==1
    figure;hold;
    plot(x,y,'ko')
    plot(min(x):max(x),feval(pfun,coefs,min(x):max(x)))
end