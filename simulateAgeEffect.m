%% Set up simulation

% These are the ages of the subjects in our sample
ages = [10.4740, 11.5362, 8.8179, 7.4928, 7.3039, 10.9888, 10.3098, ...
    7.9966, 9.2614, 9.3680, 11.1504, 12.7053, 8.6235, 7.4873, 11.5335, ...
    7.8214, 11.8594, 7.5505, 9.3299, 7.1698, 9.3680, 7.3586, 12.2672, ...
    7.3930, 7.9892, 7.4095, 7.9407, 9.9612, 10.2541, 10.2158, 12.3400];

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
simparams = nan(nrep,4,size(params,1));
r = nan(nrep,size(params,1));
for ss = 1:size(params,1)
    fprintf('\nRunning %d iterations of simulation %d\n', nrep, ss)
    % Generate a simulation of the defined effect + noise
    simdata = repmat(evalgaussian1d(params(ss,:),ages),nrep, 1);
    simnoise = randn(size(simdata)) .* noiseSD;
    simdata = simdata + simnoise;
    
    % Fit the sensitive period model to each instance of the data
    for ii = 1:nrep
        simparams(ii,:,ss) = fitgaussian1d_sd(ages, simdata(ii,:), params(ss,:));
        % Also compute correlation with age
        r(ii,ss) = corr(ages', simdata(ii,:)');
    end
end

%% Plot results

plotparams = [1 2 3 4 5 6]; % numbers to show on plots
% Calculate 68%CI for params
prc = prctile(simparams, [16 84],1);
% Extract just the SD param
prc = squeeze(prc(:,2,:));

% Plot simulated sensitive period
figure;
subplot(1,4,1); hold
x0 = 7:.1:20;
c = [linspace(0,1,size(params,1))' repmat(0,size(params,1),2)];
plotCI = [2]

for ii = 1:size(params,1)
    if any(params(ii,2) == plotCI)
        % Draw CI for a few of the simulations
        CI_x = [x0, fliplr(x0)];
        CI_y = [evalgaussian1d([params(ii,1), prc(1,ii), 1, 0],x0), ...
            fliplr(evalgaussian1d([params(ii,1), prc(2,ii), 1, 0],x0))];
        patch(CI_x, CI_y, c(ii,:),'facealpha',.4, 'edgealpha',0);
    end
    
    if any(params(ii,2) == plotparams)
        plot(x0, evalgaussian1d([params(ii,1:2) 1 0],x0),'-','color',c(ii,:),'linewidth',3);
    end
end

% Draw the period over which we have data
patch([min(ages) max(ages) max(ages) min(ages) min(ages)],[0 0 1 1 0].*1.05,[.5 .5 .5],...
    'edgealpha',.5,'facealpha',.1);

% Plot line at 2Sigma
plot([min(x0) max(x0)], repmat(evalgaussian1d([params(end,1:2) 1 0],params(end,1)+2*params(end,2)),[1 2]), '--k');

axis('tight')
axis('square')
xlabel('Age','fontsize',18); ylabel('Plastiity','fontsize',18);
set(gca, 'xtick', 7:2:20,'fontsize',14);

% Plot error on estimated parameters
subplot(1,4,2); hold
ii = 1
plot(params(ii,2), diff(prc(:,ii)),'o', 'color', c(ii,:), 'markerfacecolor', c(ii,:), 'markersize',10);

for ii = 2:size(params,1)
    % Plot unscaled as solid line
    plot([params(ii-1,2) params(ii,2)],[diff(prc(:,ii-1)) diff(prc(:,ii))],'-', 'color', c(ii,:),'linewidth',3);
    if any(params(ii,2) == plotparams)
        plot(params(ii,2), diff(prc(:,ii)),'o', 'color', 'k', 'markerfacecolor', c(ii,:), 'markersize',10);
    end
end
axis('tight')
axis('square')
set(gca, 'xtick',1:1:6, 'ygrid', 'on','fontsize',14);
xlabel('Sensitive period sigma','fontsize',18); ylabel('Estimation error','fontsize',18);

% Plot 68%CI on estimated parameters
subplot(1,4,3); hold
ii=1
patch([params(ii,2)-.2, params(ii,2)+.2, params(ii,2)+.2,params(ii,2)-.2,params(ii,2)-.2]...
    ,[prc(1,ii),prc(1,ii),prc(2,ii),prc(2,ii),prc(1,ii)]...
    ,c(ii,:),'edgecolor','k');

for ii = 2:size(params,1)
    plot([params(ii-1,2) params(ii,2)],[prc(1,ii-1) prc(1,ii)],'-', 'color', c(ii,:),'linewidth',3);
    plot([params(ii-1,2) params(ii,2)],[prc(2,ii-1) prc(2,ii)],'-', 'color', c(ii,:),'linewidth',3);
    if any(params(ii,2) == plotparams)
        patch([params(ii,2)-.2, params(ii,2)+.2, params(ii,2)+.2,params(ii,2)-.2,params(ii,2)-.2]...
            ,[prc(1,ii),prc(1,ii),prc(2,ii),prc(2,ii),prc(1,ii)]...
            ,c(ii,:),'edgecolor','k');
    end
end
axis('tight')
axis('square')

set(gca, 'xtick',1:1:6, 'ygrid', 'on','fontsize',14);
xlabel('Sensitive period sigma','fontsize',18); ylabel('68% CI on parameter estimate','fontsize',18);

% Probability of detecting a negative correlation
subplot(1,4,4); hold
rp = sum(r<0)./size(r,1);

ii = 1;
plot(params(ii,2), rp(:,1),'o', 'color', 'k', 'markerfacecolor', c(ii,:), 'markersize',10);

for ii = 2:size(params,1)
    % Plot unscaled as solid line
    plot([params(ii-1,2) params(ii,2)],[rp(:,ii-1) rp(:,ii)],'-', 'color', c(ii,:),'linewidth',3);
    if any(params(ii,2) == plotparams)
        plot(params(ii,2), rp(:,ii),'o', 'color', 'k', 'markerfacecolor', c(ii,:), 'markersize',10);
    end
end
xlabel('Sensitive period sigma','fontsize',18); ylabel('Power (Pearson r)','fontsize',18);
set(gca, 'xtick',1:1:6, 'ygrid', 'on','fontsize',14);
axis('square')

% Set size for figure
set(gcf,'position',[529 778 1630 369])
print('SensitivePeriod_sim.eps','-depsc')
print('SensitivePeriod_sim.pdf','-dpdf')

