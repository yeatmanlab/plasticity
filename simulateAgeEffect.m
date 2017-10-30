%% Set up simulation

% These are the ages of the subjects in our sample
ages = 7:.24:13;

% The noise SD is calculated based on the control subjects. Since the
% controls are not showing change over the 8 weeks, then we can calculate
% the standard deviation of change scores in the control subjects as an
% estimate of noise.
noiseSD = 0.0084;

% Mean plasticity over the age range. This is defined based on the average
% growth in the intervention subjects
mp = .0054;

% Simulation parameters
params = horzcat(repmat(min(ages), [12, 1]), [.5:.5:6]', repmat(mp, [12, 1]), zeros(12,1));
nrep = 1000; % number of iterations of simulated data

% Adjust the mean plasticity to instead reflect the max plasticity. This
% adjustment determines the peak of the gaussian.
% Justification: Huber et al (2017) measured the mean plasticity over the
% ages defined above. If we are hypothesizing that this mean came from data
% collected over a period with gaussian falloff, then this gaussian would
% have peaked much above the mean (for the younger subjects).
for ss = 1:size(params,1)
    % Calculate how much higher the peak is than the mean over this age
    % range, and then scale gaussian peak parameter appropriately
    params(ss,3) =  params(ss,3).*mp/mean(evalgaussian1d(params(ss,:),ages));
    % Test to make sure the new mean is equal to what it should be within
    % numerical precision
    assert(mean(evalgaussian1d(params(ss,:),ages)) - mp < 10^-16);
end

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
patch([min(ages) max(ages) max(ages) min(ages) min(ages)],[0 0 max(params(:,3)) max(params(:,3)) 0].*1.05,[.5 .5 .5],...
    'edgealpha',.3,'facealpha',.3);
plot([min(ages) max(ages)], [mp mp], '--k');
axis('tight')
xlabel('Age'); ylabel('Plastiity');

% Add the origin to the plot. We can assume there would be zero error for
% and SD=0;
params = vertcat([min(ages), 0, mp, 0], params);
prc = horzcat([0; 0], prc);

% Plot error on estimated parameters
subplot(1,2,2); hold
for ii = 1:size(params,1)-1
    plot([params(ii,2) params(ii+1,2)],[diff(prc(:,ii)) diff(prc(:,ii+1))],'-', 'color', c(ii,:));
    plot(params(ii+1,2), diff(prc(:,ii+1)),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
end

% Add the origin to the plot. This is
axis tight
set(gca, 'xtick',0:2:6);
grid('on')
xlabel('Sensitive period width'); ylabel('Estimation error');


return


