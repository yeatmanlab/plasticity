%% Load in data
addpath(genpath('~/git/plasticity'))
load('~/git/plasticity/data/afqOut_20180621_concat_IntSubs.mat');
rmsubs = afq.metadata.outliers |  afq.metadata.motion>0.7 ...
    | afq.sub_group==0 | afq.metadata.session>4;
afq = AFQ_RemoveSubjects(afq,rmsubs);
afq = AFQ_SubjectAvgMetadata(afq);

% Load fibers for renderings
fg = fgRead('~/git/plasticity/data/exampleFibers.mat');

%% Organize data

fgnames = AFQ_get(afq,'fgnames');
params = {'dki_MD_noden'};
nodes = 31:70;
d = table;
d.sub = afq.sub_names;
d.int_time = afq.metadata.time;
d.int_time_z = zscore(afq.metadata.time);
d.age_all = afq.metadata.visit_age;
d.age = afq.metadata.sm_visit_age;
d.age_z = zscore(afq.metadata.sm_visit_age);
d.towre_m = zscore(afq.metadata.sm_twre_index);
d.wj_m    = zscore(afq.metadata.sm_wj_brs);
d.towre_dm = afq.metadata.sdm_twre_index;
d.wj_dm = afq.metadata.sdm_wj_brs;
d.wj = afq.metadata.wj_brs;
d.towre = afq.metadata.twre_index;
d.sess = categorical(afq.metadata.session);
d.sessN = afq.metadata.session;

% Add timepoint 1 age
usubs = unique(d.sub);
for ii = 1:length(usubs)
    idx = strcmp(usubs{ii},d.sub);
    age1 = min(d.age_all(idx));
    d.age1(idx)=age1;
end
% Define young versus old subjects.
d.young = d.age1<9*12;
st = table;
stn = table;
stn_main = table;
fgnums = [1:6 9:20]

%% PCA
nc = 5; % number of pcs
md = []; fa = [];

% Collect diffusion properties in a matrix
for ii = fgnums
    md = horzcat(md,AFQ_get(afq,fgnames{ii},'dki_MD_noden'));
    %fa = horzcat(fa,AFQ_get(afq,fgnames{ii},'dki_FA_noden'));
end

% Compute PCA
[coeff, score, latent, tsquared, explained] = pca(md,'NumComponents',nc);
fprintf('\n Variance PC1=%.2f, PC2=%.2f, PC3=%.2f, PC4=%.2f, PC5=%.2f, total=%.2f',explained(1:nc),sum(explained(1:nc)))

% Add PCs into data table and fit LME
for ii = 1:nc
    d.(sprintf('pc%d',ii))=score(:,ii);
    lme = fitlme(d,sprintf('pc%d ~ int_time  + (1|sub)',ii))
    pvalmat(:,ii) = lme.Coefficients.pValue;
    lme = fitlme(d,sprintf('pc%d ~ int_time_z*age_z  + (1|sub)',ii))
    pvalmat2(:,ii) = lme.Coefficients.pValue;
end

%% Render PCA coefficients on tracts

cax = [-.03 .03]; % color range
cmap = [linspace(.1,1,128)',linspace(.1,1,128)',linspace(.8,1,128)';...
    linspace(1,.8,128)',linspace(1,.1,128)',linspace(1,.1,128)']
numf = 200;
% Loop over PCs
for pp = 1:nc
    
    % Loop over fiber groups
    fgc=0;
    for ff = fgnums
        fgc=fgc+1;
        % These are the rows containing the coeffs for fg(ff)
        rows = coeff((fgc-1)*100+1:fgc*100,:);
        
        % create color values for pc
        c = vals2colormap(rows(:,pp),cmap,cax);
        
        % Resample fibers to 100 nodes
        fgt = dtiResampleFiberGroup(fg(ff),100);
        for fff = 1:length(fgt.fibers)
            cc{fff} = c;
        end
        
        % Render
        if ff>1, nf=0; AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
        else, nf=1; lh = AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
        end
        
        clear cc c;
    end
    
    % Format and save figure
    h = colorbar;colormap(cmap);caxis(cax)
    ylabel(h,'-log_1_0(p-value) Age X Time interaction');
    axis image; axis off;
    print(sprintf('PC%dRend_1.png',pp),'-dpng','-r300');
    view(90,0); camlight(lh,'right');print(sprintf('PC%dRend_2.png',pp),'-dpng','-r300');
    view(0,90); camlight(lh,'right');print(sprintf('PC%dRend_3.png',pp),'-dpng','-r300');
    figure;h = colorbar;colormap(cmap);caxis(cax);print('cbar.eps','-depsc');
end


%% Model intervention hours x age interaction for all tract means

for ii = fgnums
    fgnospace{ii} = fgnames{ii};
    fgnospace{ii}(isspace(fgnospace{ii})) = [];
    for pp = 1:length(params)
        tmp = AFQ_get(afq,fgnames{ii},params{pp});
        d.(fgnospace{ii}) = nanmean(tmp(:,nodes),2);
        lme = fitlme(d,sprintf('%s ~ int_time*age_z  + (1|sub)',fgnospace{ii}));
        st.([fgnospace{ii} '_mean_Tstat']) = lme.Coefficients.tStat;
        st.([fgnospace{ii} '_mean_Pval']) = lme.Coefficients.pValue;
    end
end
st.Row = lme.CoefficientNames;

%% Model individual change as a function of hours

% Get the names of all the subjects
usubs = unique(d.sub);

% Create a temporary dataframe
d2 = d;
mintimes = 3; % minimum number of time points
drows = logical(zeros(size(d,1),1)); % note which rows of d were kept
% Loop over subjects
for ii = 1:length(usubs)
    % Check number of timepoints for subject ii
    if sum(strcmp(d2.sub, usubs{ii}))<3
        kr(ii) = false;
        fprintf('\nonly %d datapoints for sub %s\n',sum(strcmp(d2.sub, usubs{ii})),usubs{ii});
        d2(strcmp(d2.sub, usubs{ii}),:) = [];
    else
        % Only keep the subj if there are enough timepoints
        kr(ii) = true;
        drows = drows | strcmp(d.sub, usubs{ii}); % which rows were kept
    end
end

% Remove rows with too few time points
usubs = usubs(kr);
ind = table('RowNames',usubs);

% Fit a line to summarize each subjects rate of diffusion change over time
for ii = 1:length(usubs)
    fprintf('\nfitting slopes for sub %d',ii);
    
    % Make one entry for subject age
    ind.age(ii) = mean(d2(strcmp(d2.sub, usubs{ii}),:).age);
    
    % Loop over fiber groups
    for ff = fgnums
        fgnospace{ff} = fgnames{ff};
        fgnospace{ff}(isspace(fgnospace{ff})) = [];
        
        % loop over nodes within a fg
        tmp = AFQ_get(afq,fgnames{ff},params{pp});
        tmp = tmp(drows,:);
        % smooth rows
        for s = 1:size(tmp,1)
            tmp(s,:) = smooth(tmp(s,:),15)';
        end
        for nn = 1:100
            pf = polyfit(d2(strcmp(d2.sub, usubs{ii}),:).int_time,tmp(strcmp(d2.sub, usubs{ii}),nn),1);
            ind_nodes{ff}(ii,nn)=pf(1);
        end
        
    end
    
    % Loop over PCs
    for ff = 1:nc
        lm =  polyfit(d2(strcmp(d2.sub, usubs{ii}),:).int_time,...
            d2(strcmp(d2.sub, usubs{ii}),:).(sprintf('pc%d',ff)),1);
        ind.(['pc' num2str(ff) '_sl'])(ii) = lm(1);
    end
end
% Rescale MD units
ind{:,2:end} = -1000.*ind{:,2:end};

% Report stats relating individual PC change and age
for ff = 1:nc
    [rval, pval] = corr(ind.age,ind.(['pc' num2str(ff) '_sl']));
    fprintf('\nPC%d, r=%.2f, p=%.2f, BF=%.2f',ff,rval,pval,corrbf(rval,length(ind.age)));
    [~,pval,~,stats] = ttest2(ind(ind.age<=108,:).(['pc' num2str(ff) '_sl']),...
        ind(ind.age>108,:).(['pc' num2str(ff) '_sl']));
    fprintf('\nPC%d, t(%d)=%.2f, p=%.2f',ff,stats.df,stats.tstat,pval);
    [~,pval,~,stats] = ttest(ind.(['pc' num2str(ff) '_sl']),0);
    fprintf('\nPC%d,T-test against zero change t(%d)=%.2f, p=%.2f\n',ff,stats.df,stats.tstat,pval);
end

%% Redo analysis just looking at first 2 timepoints

% Get the names of all the subjects
usubs = unique(d.sub);

% Create another temporary dataframe
d3 = d;
drows = logical(zeros(size(d,1),1)); % note which rows of d were kept
% Loop over subjects
for ii = 1:length(usubs)
    % Check to make sure the subject has data for session 1 and 2
    if ~any(d3(strcmp(d3.sub, usubs{ii}),:).sessN == 1) || ~any(d3(strcmp(d3.sub, usubs{ii}),:).sessN == 2)
        kr(ii) = false;
        fprintf('\nFor sub %s sessions:',usubs{ii});
        [d3(strcmp(d3.sub, usubs{ii}),:).sessN']
        d3(strcmp(d3.sub, usubs{ii}),:) = [];
    else
        % Only keep the subj if there are enough timepoints
        kr(ii) = true;
        drows = drows | strcmp(d.sub, usubs{ii}); % which rows were kept
    end
end

% Remove rows with too few time points
usubs2 = usubs(kr);
ind2 = table('RowNames',usubs2);
for ii = 1:length(usubs2)
    
    % add initial age
    ind2.age(ii) = d3(strcmp(d3.sub, usubs2{ii}),:).age1(1);
    % Loop over PCs
    for ff = 1:nc
        lm =  polyfit(d3(strcmp(d3.sub, usubs2{ii}),:).int_time,...
            d3(strcmp(d3.sub, usubs2{ii}),:).(sprintf('pc%d',ff)),1);
        ind2.(['pc' num2str(ff) '_sl'])(ii) = lm(1);
        % Fit model for 2 sessions
        lm2 =  polyfit(d3(strcmp(d3.sub, usubs2{ii}),:).int_time(1:2),...
            d3(strcmp(d3.sub, usubs2{ii}),:).(sprintf('pc%d',ff))(1:2),1);
        ind2.(['pc' num2str(ff) '_sl2'])(ii) = lm2(1);
        
        % Compute percentage of change between 60 and 160 hours. We compute
        % the overall change between 0 and 160 hours based on the model
        % taking into account 4 data points. We compute the change between
        % 0 and 60 hours based on the 2 timepoint model. The reason for
        % using the model rather than the subtraction is to control for
        % differences in the timing of sessions among subjects.
        ind2.(['pc' num2str(ff) 'perc_1_2'])(ii) = ...
            (polyval(lm,160) - polyval(lm2,60)) ./ (polyval(lm,160) - polyval(lm,0));
    end
end

% Compute correlation between age and pc change
for pcn = 1:5
    [rval, pval] = corr(ind2.age,ind2.(sprintf('pc%dperc_1_2',pcn)),'type','spearman');
    fprintf('\nPC%d session 1-2 percent, r = %.2f, p = %.2f',pcn, rval, pval);
    [rval, pval] = corr(ind2.age,ind2.(sprintf('pc%d_sl2',pcn)),'type','spearman');
    fprintf('\nPC%d session 1-2, r = %.2f, p = %.2f',pcn, rval, pval);
end
%% Model change over intervention hours x age interaction for each node
dn = d; % New table to contain each node
for ii = fgnums
    fgnospace{ii} = fgnames{ii};
    fgnospace{ii}(isspace(fgnospace{ii})) = [];
    fprintf('\n Modeling plasticity in %s\n',fgnospace{ii});
    for pp = 1:length(params)
        tmp = AFQ_get(afq,fgnames{ii},params{pp});
        % smooth rowas
        for s = 1:size(tmp,1)
            tmp(s,:) = smooth(tmp(s,:),15)';
        end
        for nn = 1:100
            % Simple model (no interaction)
            dn.([fgnospace{ii} '_node' num2str(nn)]) = tmp(:,nn);
            lme = fitlme(dn,sprintf('%s ~ int_time + (1|sub)',[fgnospace{ii} '_node' num2str(nn)]));
            stn_main.([fgnospace{ii} '_node' num2str(nn) '_Tstat']) = lme.Coefficients.tStat;
            stn_main.([fgnospace{ii} '_node' num2str(nn) '_Pval']) = lme.Coefficients.pValue;
            % Model with interaction
            lme = fitlme(dn,sprintf('%s ~ int_time*age_z  + (1|sub)',[fgnospace{ii} '_node' num2str(nn)]));
            stn.([fgnospace{ii} '_node' num2str(nn) '_Tstat']) = lme.Coefficients.tStat;
            stn.([fgnospace{ii} '_node' num2str(nn) '_Pval']) = lme.Coefficients.pValue;
        end
    end
end
stn.Row = lme.CoefficientNames;

%% Make rendering
numf = 200; % number of fibers to render

% color fibers based on correlation (pearson's r) between age and MD change
cax = [-.8 .8]; % color range
cmap =AFQ_colormap('bgr',256);
clear rval pval;
for ff = fgnums
    [~,t_pval,~,stats] = ttest(ind_nodes{ff});
    [rval(ff,:), pval] = corr(ind.age, ind_nodes{ff});
    fprintf('\n%s max r = %.2f, p = %.3f\n',fgnames{ff},max(rval(ff,:)),min(pval))
    c = vals2colormap(rval(ff,:),cmap,cax);
    fgt = dtiResampleFiberGroup(fg(ff),100);
    for fff = 1:length(fgt.fibers)
        cc{fff} = c;
    end
    
    if ff>1, nf=0; AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
    else, nf=1; lh = AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
    end
    
    clear cc c;
end
h = colorbar;colormap(cmap);caxis(cax)
ylabel(h,'Correlation with age');
axis image; axis off;
print('Rvals1.png','-dpng','-r300');
view(90,0); camlight(lh,'infinite');print('Rvals2.png','-dpng','-r300');
view(0,90); camlight(lh,'infinite');print('Rvals3.png','-dpng','-r300');
figure;h = colorbar;colormap(cmap);caxis(cax);print('cbar_Rvals.eps','-depsc');
% Make histogram figure
figure;fgn=0;
for ff = fgnums
    fgn = fgn+1;
    subplot(3,6,fgn);histogram(rval(ff,:),20);
    title(fgnames{ff});xlabel('r value');
end
set(gcf,'position',[0 100 1000 600]);
print('RvalueHist.eps','-depsc');

% Color fibers based on p-value for Time X Age interaction (from LME)
cax = [0 3]; % color range
cmap = [linspace(.3,1,255)',linspace(.3,1,255)',ones(255,1).*.3];
fgc=0;
for ff = fgnums
    fgc=fgc+1;
    % These are the columns containing the p-value for fg(ff)
    cols = (fgc-1)*200+2:2:fgc*200;
    pval(ff,:) = table2array(stn(4,cols));
    fprintf('\n%s AgeXTime Interaction min p = %.3f\n',fgnames{ff},min(pval(ff,:)))
    
    % Convert table to array and -log10 it
    logp(ff,:) = -log10(pval(ff,:));
    
    c = vals2colormap(logp(ff,:),cmap,cax);
    fgt = dtiResampleFiberGroup(fg(ff),100);
    for fff = 1:length(fgt.fibers)
        cc{fff} = c;
    end
    
    if ff>1, nf=0; AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
    else, nf=1; lh = AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
    end
    
    clear cc c;
end
h = colorbar;colormap(cmap);caxis(cax)
ylabel(h,'-log_1_0(p-value) Age X Time interaction');
axis image; axis off;
print('Pvals_Int1.png','-dpng','-r300');
view(90,0); camlight(lh,'infinite');print('Pvals_Int2.png','-dpng','-r300');
view(0,90); camlight(lh,'infinite');print('Pvals_Int3.png','-dpng','-r300');
figure;h = colorbar;colormap(cmap);caxis(cax);print('cbar_Int.eps','-depsc');
% Make histogram figure
figure;fgn=0;
for ff = fgnums
    fgn = fgn+1;
    subplot(3,6,fgn);h=histogram(logp(ff,:),20,'BinLimits',[0 4]);axis('tight')
    text(h.BinEdges(max(find(h.Values)))+h.BinWidth,max(h.Values)/15,...
        ['\downarrow' sprintf('p = %.3f',min(pval(ff,:)))],'color',[.7 0 0]);
    title(fgnames{ff});xlabel('p-value');
    set(gca,'xtick',[2 4],'xticklabels',{'0.01' '0.0001'});
end
set(gcf,'position',[0 100 1000 600]);
print('PvalueHist_Interaction.eps','-depsc');

% Color fibers based on p-value for main effect of Time (from LME)
cax = [0 3]; % color range
fgc=0;
for ff = fgnums
    fgc=fgc+1;
    % These are the columns containing the p-value for fg(ff)
    cols = (fgc-1)*200+2:2:fgc*200;
    pval(ff,:) = table2array(stn_main(2,cols));
    fprintf('\n%s Main Effect Time min p = %.3f\n',fgnames{ff},min(pval(ff,:)))
    
    % Convert table to array and -log10 it
    logp(ff,:) = -log10(pval(ff,:));
    
    c = vals2colormap(logp(ff,:),cmap,cax);
    fgt = dtiResampleFiberGroup(fg(ff),100);
    for fff = 1:length(fgt.fibers)
        cc{fff} = c;
    end
    
    if ff>1, nf=0; AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
    else, nf=1; lh = AFQ_RenderFibers(fgt,'numfibers',numf,'color',cc,'newfig',nf);
    end
    
    clear cc c;
end
h = colorbar;colormap(cmap);caxis(cax)
ylabel(h,'-log_1_0(p-value) Change over intervention');
axis image; axis off;
print('Pvals_time1.png','-dpng','-r300');
view(90,0); camlight(lh,'infinite');print('Pvals_time2.png','-dpng','-r300');
view(0,90); camlight(lh,'infinite');print('Pvals_time3.png','-dpng','-r300');
% Make histogram figure
figure;fgn=0;
for ff = fgnums
    fgn = fgn+1;
    subplot(3,6,fgn);h=histogram(logp(ff,:),20,'BinLimits',[0 4]);axis('tight')
    text(h.BinEdges(max(find(h.Values)))+h.BinWidth,max(h.Values)/15,...
        [sprintf('p = %.4f',min(pval(ff,:))) '\downarrow'],'color',[.7 0 0], 'horizontalAlignment', 'right');
    title(fgnames{ff});xlabel('p-value');
    set(gca,'xtick',[2 4],'xticklabels',{'0.01' '0.0001'});
end
set(gcf,'position',[0 100 1000 600]);
print('PvalueHist_MainEffect.eps','-depsc');

%% Show growth trajectories for younger versus older subjects

% Fit separate models for PC1
lme = fitlme(d,'pc1 ~ sess + (1 | sub)')
lmeY = fitlme(d(d.young,:),'pc1 ~ sess + (1 | sub)');
lmeO = fitlme(d(~d.young,:),'pc1 ~ sess + (1 | sub)');

% Plot results
figure;
subplot(1,2,1);hold
errorbar([lme.Coefficients.Estimate(1), lme.Coefficients.Estimate(2:end)'+lme.Coefficients.Estimate(1)],...
    lme.Coefficients.SE, '-o','linewidth',3,'color', [.5 0 .5],'markerfacecolor', [.5 0 .5]);
plot([lmeY.Coefficients.Estimate(1), lmeY.Coefficients.Estimate(2:end)'+lmeY.Coefficients.Estimate(1)],...
    '-o','linewidth',3,'color', [.7 0 0],'markerfacecolor', [.7 0 0]);
plot([lmeO.Coefficients.Estimate(1), lmeO.Coefficients.Estimate(2:end)'+lmeO.Coefficients.Estimate(1)],...
    '-o','linewidth',3,'color', [0 0 .7],'markerfacecolor', [0 0 .7]);
grid on;set(gca,'gridalpha',.4)
set(gca,'fontsize',14)
xlabel('Session','fontname','Helvetica','fontsize',16),ylabel('Mean diffusivity (\mum^2/ms)','fontname','Helvetica','fontsize',16)


% Repeate for behavior
lme = fitlme(d,'wj ~ sess + (1 | sub)')
lmeY = fitlme(d(d.young,:),'wj ~ sess + (1 | sub)');
lmeO = fitlme(d(~d.young,:),'wj ~ sess + (1 | sub)');

% Plot results
subplot(1,2,2);hold
h(1) = errorbar([lme.Coefficients.Estimate(1), lme.Coefficients.Estimate(2:end)'+lme.Coefficients.Estimate(1)],...
    lme.Coefficients.SE, '-o','linewidth',3,'color', [.5 0 .5],'markerfacecolor', [.5 0 .5]);
h(2) = plot([lmeY.Coefficients.Estimate(1), lmeY.Coefficients.Estimate(2:end)'+lmeY.Coefficients.Estimate(1)],...
    '-o','linewidth',3,'color', [.7 0 0],'markerfacecolor', [.7 0 0]);
h(3) = plot([lmeO.Coefficients.Estimate(1), lmeO.Coefficients.Estimate(2:end)'+lmeO.Coefficients.Estimate(1)],...
    '-o','linewidth',3,'color', [0 0 .7],'markerfacecolor', [0 0 .7]);
grid on;set(gca,'gridalpha',.4)
set(gca,'fontsize',14)
xlabel('Session','fontname','Helvetica','fontsize',16);ylabel('Basic Reading Skills (standard score)','fontname','Helvetica','fontsize',16);
legend(h,{'All' 'Younger' 'Older'},'location','NorthWest')
set(gcf, 'Position',  [100, 100, 600, 350]);

% Save
print('Figure3_GrowthCurves.eps','-depsc')

%% Does behavioral growth depend on age?
lme_wj = fitlme(d,'wj ~ int_time_z * age_z + (1|sub)')
lme_tw = fitlme(d,'towre ~ int_time_z * age_z + (1|sub)')

return

%% Plots of plasticity versus age
figure; hold
age = ind.age./12;
lme = fitlme(d,'pc1 ~ int_time + (int_time | sub)');
[B,Bnames,RE] = randomEffects(lme)
RE(~strcmp(stats.Name,'int_time'),:)= [];
% This is an estimnate of each individuals growth rate. NOTE WE MAKE IT
% NEGATIVE FOR PLOTTING PURPOSES
m_est = lme.Coefficients.Estimate(2); % average growth rate
m_est_se = lme.Coefficients.SE(2); % SE on average growth rate
ind_est = -(RE.Estimate + m_est);
% Get the age for each subject
for ii = 1:size(RE,1)
    % Get age (scaled to years /12) for each sub
    lme_age(ii) = d(strcmp(d.sub,RE.Level{ii}),:).age(1)./12;
end

% baseed on individual estimate
plot(age, ind.pc1_sl,'ko','markerfacecolor','k')
lsline

% baseed on estimates from LME
figure; hold
plot(lme_age, ind_est,'ko','markerfacecolor','k')
lsline

% Plot gaussian model
xx = 6:.1:13;
plot(xx,evalgaussian1d([6 1 max(ind_est).*.8 -m_est],xx))

%% Fit gaussian model to data

for ff = fgnums
    lm =  fitlm(ind,sprintf('%s_node1_sl ~ age^2',fgnospace{ff}));
    lm_R2(ff) = lm.Rsquared.Ordinary.*100;
    % Seed params for gaussian fitting
    [x1,x2] = meshgrid(1:10, 0:.001: max(ind.(sprintf('%s_node1_sl',fgnospace{ff})))); params0 = [x1(:), x2(:)];
    params0 = [repmat(min(ind.age),size(params0,1),1), params0, zeros(size(params0,1),1)];
    % fit the gaussian
    [params(ff,:), R2(ff), params_tmp, R2_tmp]= fitgaussian1d_sd(ind.age,ind.(sprintf('%s_node1_sl',fgnospace{ff})),...
        params0,'sgd');
    
    figure; hold('on')
    plot(ind.age, ind.(sprintf('%s_node1_sl',fgnospace{ff})),'o',...
        'color',[.3 .6 .8],'markerfacecolor',[.3 .6 .8]);
    axis tight
    L(1) = lsline;set(L(1),'linewidth',2,'color',[.3 .6 .8]);
    L(2) = refline(0,0);set(L(2),'linewidth',1,'color',[0 0 0],'linestyle','--');
    text(mean(ind.age),prctile(ind.(sprintf('%s_node1_sl',fgnospace{ff})),90),...
        sprintf('T-int = %.2f, T-slope = %.2f',lm.Coefficients.tStat))
    plot(linspace(min(ind.age),max(ind.age),500),evalgaussian1d(params(ff,:),linspace(min(ind.age),max(ind.age),500)),'-r')
    hold('off')
end
