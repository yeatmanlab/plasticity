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
params = horzcat(repmat(min(ages), [26, 1]), [1:.2:6]', repmat(mp, [26, 1]), zeros(26,1));
nrep = 1000; % number of iterations of simulated data

% Adjust the mean plasticity to instead reflect the max plasticity. This
% adjustment determines the peak of the gaussian.
% Justification: Huber et al (2017) measured the mean plasticity over the
% ages defined above. If we are hypothesizing that this mean came from data
% collected over a period with gaussian falloff, then this gaussian would
% have peaked much above the mean (for the younger subjects).
scaleGaussian = true;
if scaleGaussian
    params = scalePeak(params, mp, ages);
end


%% Run simulation

for ss = 1:size(params,1)
    fprintf('\nRunning %d iterations of simulation %d\n', nrep, ss)
    % Generate a simulation of the defined effect + noise
    simdata = repmat(evalgaussian1d(params(ss,:),ages),nrep, 1);
    simnoise = randn(size(simdata)) .* noiseSD;
    simdata = simdata + simnoise;
    simdata_scaled = repmat(evalgaussian1d(params_scaled(ss,:),ages),nrep, 1);
    % Generate independent noise
    simnoise = randn(size(simdata_scaled)) .* noiseSD;
    simdata_scaled = simdata_scaled + simnoise;
    
    % Fit the sensitive period model to each instance of the data
    for ii = 1:nrep
        simparams(ii,:,ss) = fitgaussian1d_sd(ages, simdata(ii,:), params(ss,:));
        simparams_scaled(ii,:,ss) = fitgaussian1d_sd(ages, simdata_scaled(ii,:), params_scaled(ss,:));
        % Also compute correlation with age
        r(ii,ss) = corr(ages', simdata(ii,:)');
        r_scaled(ii,ss) = corr(ages', simdata_scaled(ii,:)');
    end
end

%% Plot results

plotparams = [1 2 3 4 5 6]; % numbers to show on plots
% Calculate 68%CI for params
prc = prctile(simparams, [16 84],1);
prc_scaled = prctile(simparams_scaled, [16 84],1);
% Extract just the SD param
prc = squeeze(prc(:,2,:));
prc_scaled = squeeze(prc_scaled(:,2,:));

% Plot simulated sensitive period
figure;
subplot(1,4,1); hold
x0 = 7:.1:20;
c = parula(size(params,1));
for ii = 1:size(params,1)
    if any(params(ii,2) == plotparams)
        plot(x0, evalgaussian1d(params(ii,:),x0),'-','color',c(ii,:),'linewidth',3);
    end
end
patch([min(ages) max(ages) max(ages) min(ages) min(ages)],[0 0 max(params(:,3)) max(params(:,3)) 0].*1.05,[.5 .5 .5],...
    'edgealpha',.3,'facealpha',.3);
plot([min(ages) max(ages)], [mp mp], '--k');
axis('tight')
xlabel('Age'); ylabel('Plastiity');

% Plot error on estimated parameters
subplot(1,4,2); hold
ii = 1
plot(params(ii,2), diff(prc(:,ii)),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
plot(params_scaled(ii,2), diff(prc_scaled(:,ii)),'o', 'color', [0 0 0], 'markerfacecolor', c(ii,:));

for ii = 2:size(params,1)
    % Plot unscaled as solid line
    plot([params(ii-1,2) params(ii,2)],[diff(prc(:,ii-1)) diff(prc(:,ii))],'-', 'color', c(ii,:),'linewidth',3);
    plot([params_scaled(ii-1,2) params_scaled(ii,2)],[diff(prc_scaled(:,ii-1)) diff(prc_scaled(:,ii))],'-', 'color', c(ii,:),'linewidth',3);
    if any(params(ii,2) == plotparams)
        plot(params(ii,2), diff(prc(:,ii)),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
        plot(params_scaled(ii,2), diff(prc_scaled(:,ii)),'o', 'color', [0 0 0], 'markerfacecolor', c(ii,:));
    end
end
axis tight
set(gca, 'xtick',0:2:6);
grid('on')
xlabel('Sensitive period width'); ylabel('Estimation error');

% Plot 68%CI on estimated parameters
subplot(1,4,3); hold

for ii = 2:size(params,1)
    plot([params_scaled(ii-1,2) params_scaled(ii,2)],[prc_scaled(1,ii-1) prc_scaled(1,ii)],'-', 'color', c(ii,:),'linewidth',3);
    plot([params_scaled(ii-1,2) params_scaled(ii,2)],[prc_scaled(2,ii-1) prc_scaled(2,ii)],'-', 'color', c(ii,:),'linewidth',3);
    if any(params(ii,2) == plotparams)
        patch([params_scaled(ii,2)-.1, params_scaled(ii,2)+.1, params_scaled(ii,2)+.1,params_scaled(ii,2)-.1,params_scaled(ii,2)-.1]...
            ,[prc_scaled(1,ii),prc_scaled(1,ii),prc_scaled(2,ii),prc_scaled(2,ii),prc_scaled(1,ii)]...
            ,c(ii,:),'edgecolor',c(ii,:));
    end
end
axis tight
set(gca, 'xtick',0:2:6);
grid('on')
xlabel('Sensitive period width'); ylabel('68% CI on parameter estimate');

% Probability of detecting a negative correlation
subplot(1,4,4); hold
rp = sum(r<0)./size(r,1);
rp_scaled = sum(r_scaled<0)./size(r_scaled,1);

ii = 1;
plot(params(ii,2), rp(:,1),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
plot(params(ii,2), rp_scaled(:,ii),'o', 'color', [0 0 0], 'markerfacecolor', c(ii,:));

for ii = 2:size(params,1)
    % Plot unscaled as solid line
    plot([params(ii-1,2) params(ii,2)],[rp(:,ii-1) rp(:,ii)],'-', 'color', c(ii,:),'linewidth',3);
    plot([params(ii-1,2) params(ii,2)],[rp_scaled(:,ii-1) rp_scaled(:,ii)],'-', 'color', c(ii,:),'linewidth',3);
    if any(params(ii,2) == plotparams)
        plot(params(ii,2), rp(:,ii),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:));
        plot(params(ii,2), rp_scaled(:,ii),'o', 'color', [0 0 0], 'markerfacecolor', c(ii,:));
    end
end
xlabel('Sensitive period width'); ylabel('Power');
grid('on')


