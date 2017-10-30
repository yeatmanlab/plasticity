%% Set up simulation

% These are the ages of the subjects in our sample
ages = 7:.24:13

% The noise SD is calculated based on the control subjects. Since the
% controls are not showing change over the 8 weeks, then we can calculate
% the standard deviation of change scores in the control subjects as an
% estimate of noise.
noiseSD = .4

% Mean plasticity over the age range. This is defined based on the average
% growth in the intervention subjects
mp = 1;

% Simulation parameters
params = horzcat(repmat(min(ages), [11, 1]), [1:.5:6]', repmat(mp, [11, 1]), zeros(11,1));
nrep = 1000; % number of iterations of simulated data

% Adjust the mean plasticity to instead reflect the max plasticity. This
% adjustment determines the peak of the gaussian.
% TODO

%% Run simulation

for ss = 1:size(params,1)
    fprintf('\nRunning %d iterations of simulation %d\n', nrep, ss)
    % Generate a simulation of the defined effect + noise
    simdata = repmat(evalgaussian1d(params(ss,:),ages),nrep, 1);
    simnoise = randn(size(simdata)) .* noiseSD;
    simdata = simdata + simnoise;
    
    % Fit the sensitive period model to each instance of the data
    for ii = 1:nrep
        simparams(ii,:,ss) = fitgaussian1d_sd(ages, simdata(ii,:), params(ss,:));
    end
end

%% Plot results

% Calculate 68%CI for params
prc = prctile(simparams, [16 84],1);
% Extract just the SD param
prc = squeeze(prc(:,2,:));

% Plot simulated sensitive period
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

% Plot error on estimated parameters
subplot(1,2,2); hold
for ii = 1:size(params,1)-1
    plot([params(ii,2) params(ii+1,2)],[diff(prc(:,ii)) diff(prc(:,ii+1))],'-', 'color', c(ii,:));
    plot(params(ii,2), diff(prc(:,ii)),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
end
plot(params(end,2), diff(prc(:,end)),'o', 'color', c(end,:), 'markerfacecolor', c(end,:));

axis([0 6 0 24]);
xlabel('Sensitive period width'); ylabel('Estimation error');


return 


