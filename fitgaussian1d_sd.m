function [params,R2] = fitgaussian1d_sd(x,y,params0)

% function [params,R2] = fitgaussian1d(x,y,params0)
%
% <x>,<y> are row vectors of the same size.  <x> specifies x-coordinates
%   and <y> specifies values.  if <x> is [], default to 1:length(<y>).
% <params0> is an initial seed. We only fit the SD so the other params get
% fixed as they are set in params0
%
% use lsqcurvefit.m to estimate parameters of a 1D Gaussian function.
% return:
%  <params> is like the input to evalgaussian1d.m
%  <R2> is the R^2 between fitted and actual y-values (see calccod.m).
%
% example:
% xx = 1:.1:10;
% yy = evalgaussian1d([5 1 4 0],xx);
% yy = yy + 0.1*randn(1,length(yy));
% [params,R2] = fitgaussian1d(xx,yy);
% figure; hold on;
% plot(xx,yy,'ro-');
% plot(xx,evalgaussian1d(params,xx),'b-');

% input

% construct coordinates
if isempty(x)
  x = 1:length(y);
end

% define options
options = optimset('Display','off','FunValCheck','on','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun',1e-6,'TolX',1e-6);

% define bounds
%              m    s    g    d
paramslb = 0;
paramsub = Inf;

% do it
[params,d,d,exitflag,output] = lsqcurvefit(@(pp,xx) ...
  evalgaussian1d([params0(1) pp params0(3:4)],xx),params0(2),x,y,paramslb,paramsub,options);
assert(exitflag > 0);
params = [params0(1) params params0(3:4)];

% how well did we do?
R2 = calccod(evalgaussian1d(params,x),y);
