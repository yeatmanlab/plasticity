%%
ages = 7:.24:13
noiseSD = .4
params = [7 1 1 0; 7 2 1 0; 7 3 1 0; 7 4 1 0; 7 5 1 0; 7 6 1 0]
nrep = 100;

for ss = 1:size(params,1)
    simdata = repmat(evalgaussian1d(params(ss,:),ages),nrep, 1);
    simnoise = randn(size(simdata)).*noiseSD;
    simdata = simdata + simnoise;
    for ii = 1:nrep
        simparams(ii,:,ss) = fitgaussian1d_sd(ages, simdata(ii,:), params(ss,:));
    end
end

%% Plot
std_simparams = squeeze(std(simparams));
figure;
subplot(1,2,1); hold
x0 = 7:.1:20;
c = parula(size(params,1));
for ii = 1:size(params,1)
    plot(x0, evalgaussian1d(params(ii,:),x0),'-','color',c(ii,:));
end
patch([min(ages) max(ages) max(ages) min(ages) min(ages)],[0 0 1 1 0],[.5 .5 .5],...
    'edgealpha',.3,'facealpha',.3)
axis('tight')
xlabel('Age'); ylabel('Plastiity');
subplot(1,2,2); hold
for ii = 1:size(params,1)-1
    plot([params(ii,2) params(ii+1,2)],[std_simparams(2,ii) std_simparams(2,ii+1)],'-', 'color', c(ii,:));
    plot(params(ii,2),std_simparams(2,ii),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
end
plot(params(end,2),std_simparams(2,end),'o', 'color', c(end,:), 'markerfacecolor', c(end,:));

axis([0 6 0 24]);
xlabel('Sensitive period width'); ylabel('Estimation error');
% 
% subplot(2,2,3); hold
% for ii = 1:size(params,1)
%     plot(x0, evalPoissonCurve(paramsPois(ii,:),x0),'-','color',c(ii,:));
% end
% subplot(2,2,2); hold
% for ii = 1:size(params,1)-1
%     plot([params(ii,2) params(ii+1,2)],[std_simparams(2,ii) std_simparams(2,ii+1)],'-', 'color', c(ii,:));
%     plot(params(ii,2),std_simparams(2,ii),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
% end

% parameter distributsion

return 


