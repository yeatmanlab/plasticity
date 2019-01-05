function [params, R2, params_tmp, R2_tmp] = fitgaussian1d_sd(x,y,params0,fitparams)

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

if ~exist('fitparams','var') || isempty(fitparams)
    fitparams = 's';
end
% define options
options = optimset('Display','off','FunValCheck','on','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun',1e-12,'TolX',1e-12);

% define bounds
%              m    s    g    d

% Loop over starting params and find best fit. This avoids local minimum
for ii = 1:size(params0,1)
    if strcmp(fitparams,'s')
        paramslb = 0;
        paramsub = Inf;
        [paramsFit,d,d,exitflag,output] = lsqcurvefit(@(pp,xx) ...
            evalgaussian1d([params0(ii,1) pp params0(ii,3:4)],xx),params0(ii,2),x,y,paramslb,paramsub,options);
        assert(exitflag > 0);
        params_tmp(ii,:) = [params0(ii,1) paramsFit params0(ii,3:4)];
    elseif strcmp(fitparams,'sg')
        paramslb = [0 0];
        paramsub = [Inf  Inf];
        [paramsFit,d,d,exitflag,output] = lsqcurvefit(@(pp,xx) ...
            evalgaussian1d([params0(ii,1) pp params0(ii,4)],xx),params0(ii,2:3),x,y,paramslb,paramsub,options);
        assert(exitflag > 0);
        params_tmp(ii,:) = [params0(ii,1) paramsFit params0(ii,4)];
    elseif strcmp(fitparams,'sgd')
        paramslb = [0 0 -Inf];
        paramsub = [Inf  Inf  Inf];
        [paramsFit,d,d,exitflag,output] = lsqcurvefit(@(pp,xx) ...
            evalgaussian1d([params0(ii,1) pp],xx),params0(ii,2:4),x,y,paramslb,paramsub,options);
        assert(exitflag > 0);
        params_tmp(ii,:) = [params0(ii,1) paramsFit];
    end
    
    % how well did we do?
    R2_tmp(ii) = calccod(evalgaussian1d(params_tmp(ii,:),x),y);
end
[~,ix]=max(R2_tmp);
params = params_tmp(ix,:);
R2 = R2_tmp(ix);