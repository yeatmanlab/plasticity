%% Load in data
addpath(genpath('~/git/plasticity'))
load('~/git/plasticity/data/afqOut_20190715_meta_motion.mat');
% load('~/git/plasticity/data/motdisprms.mat');

% See if there are any subjects with crazy outlier values
afq.metadata.outliers = AFQ_outliers(afq, {'dki_MD_nnk'}, 4, 40);

% Remove subjects with lots of motion or lots of outliers
rmsubs = afq.metadata.outliers |  afq.metadata.motion>3 ...
    | afq.metadata.session>4
afq = AFQ_RemoveSubjects(afq,rmsubs);
afq = AFQ_SubjectAvgMetadata(afq);


%% Organize data

d = table;
d.sub = afq.sub_names;
d.int_days = afq.metadata.int_time-nanmean(afq.metadata.int_time); % center hours
d.int_time = afq.metadata.int_hours-nanmean(afq.metadata.int_hours); % center hours
d.int_time_z = zscore(afq.metadata.int_hours);
d.age_all = afq.metadata.visit_age;
d.age = afq.metadata.sm_visit_age-nanmean(afq.metadata.sm_visit_age); % center age
d.age_z = zscore(afq.metadata.sm_visit_age);
d.towre_m = zscore(afq.metadata.sm_twre_index);
d.wj_m    = zscore(afq.metadata.sm_wj_brs);
d.towre_dm = afq.metadata.sdm_twre_index;
d.wj_dm = afq.metadata.sdm_wj_brs;
d.wj = afq.metadata.wj_brs;
d.towre = afq.metadata.twre_index;
d.sess = categorical(afq.metadata.session);
d.sessN = afq.metadata.session;
d.sub_group = categorical(afq.sub_group);

% Make a table to save results
st = table;
st1=table;
st0=table;
fgnums = [1:6 9:20];

fgnames = AFQ_get(afq,'fgnames');
params = {'dki_MD_nnk'};
nodes = 21:80;
% abbreviated names for plot labels
tractnames = {'LThal','RThal','LCST','RCST','LCingC','RCingC','LCingH','RCingH',...
    'CFMajor','CFMinor','LIFOF','RIFOF','LILF','RILF','LSLF','RSLF','LUnc','RUnc','LArc','RArc'};


%% TIME X GROUP INTERACTION
figure; hold
for ii = fgnums
    fgnospace{ii} = fgnames{ii};
    fgnospace{ii}(isspace(fgnospace{ii})) = [];
    for pp = 1:length(params)
        tmp = AFQ_get(afq,fgnames{ii},params{pp});
        d.(fgnospace{ii}) = nanmean(tmp(:,nodes),2);
        lme0 = fitlme(d(d.sub_group=='0',:),sprintf('%s ~ int_days + (1|sub)',fgnospace{ii}),'dummyvarcoding','effect');
        lme1 = fitlme(d(d.sub_group=='1',:),sprintf('%s ~ int_days + (1|sub)',fgnospace{ii}),'dummyvarcoding','effect');
        lme = fitlme(d,sprintf('%s ~ int_days * sub_group + (1|sub)',fgnospace{ii}),'dummyvarcoding','effect');
        st.([fgnospace{ii} '_mean_Tstat']) = lme.Coefficients.pValue;
        st0.([fgnospace{ii} '_mean_Pval']) = lme0.Coefficients.pValue;
        st1.([fgnospace{ii} '_mean_Pval']) = lme1.Coefficients.pValue;
        h(1) = errorbar(ii,lme0.Coefficients.Estimate(2),lme0.Coefficients.SE(2),'ko');
        h(2) =  errorbar(ii,lme1.Coefficients.Estimate(2),lme1.Coefficients.SE(2),'ro');
        
    end
end
legend(h,{'control' 'intervention'}); grid on
st.Row = lme.CoefficientNames;
