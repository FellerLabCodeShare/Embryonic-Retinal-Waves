% Created by Christiane 9/29/21
%% This code 
% 1) loads the inter-event-interval table (table with % of ROIs participating in each wave. since ROIs tiled the whole retina we used this percentage as a proxy for measuring wave size. 
% 2) Separate data by condition
% 3) Run the stats
% 4) plots a histogram with the mean value

%% 1) Load neuronTable
iwiVariable_path = 'C:\Users\Christiane\Dropbox\Lab_Data\ipRGC_participation_in_waves\4. CompiledResults\Stage 1 waves\SpatioTemp\Matlab variables for iwi';
cd(iwiVariable_path);
load('neuronIwiTable.mat');

%% 2) normalize iwis by experiment and condition
edges = 0:15:200;
% wt ctrl
numCells_table = unique(neuronTable.totalNeuronNum,'stable');
ctrl_iwi = neuronTable.ctrl_meanWaveIntervals;
ctrl_iwi(isnan(ctrl_iwi))=[];
r_1 = neuronTable(strcmp(neuronTable.movie_name, '211119_n3_r4_concatenated stack_motionCorr_dfof'),:);
r1 = r_1.ctrl_meanWaveIntervals;
r1(isnan(r1))=[];

r_2 = neuronTable(strcmp(neuronTable.movie_name, '211119_n4_r6_concatenated stack_motionCorr_dfof'),:);
r2 = r_2.ctrl_meanWaveIntervals;
r2(isnan(r2))=[];

r_3 = neuronTable(strcmp(neuronTable.movie_name, '211119_n6_r9_concatenated stack_motionCorr_dfof'),:);
r3 = r_3.ctrl_meanWaveIntervals;
r3(isnan(r3))=[];

r_4 = neuronTable(strcmp(neuronTable.movie_name, '211119_n6_r10_concatenated stack_motionCorr_dfof'),:);
r4 = r_4.ctrl_meanWaveIntervals;
r4(isnan(r4))=[];

% b2 ctrl
ctrlB2_iwi = neuronTable.ctrlB2_meanWaveIntervals;
ctrlB2_iwi(isnan(ctrlB2_iwi))=[];

b_1 = neuronTable(strcmp(neuronTable.movie_name, '220215_n2_r1_concatenated_motionCorr_dfof'),:);
b1 = b_1.ctrlB2_meanWaveIntervals;
b1(isnan(b1))=[];

b_2 = neuronTable(strcmp(neuronTable.movie_name, '220215_n4_r1_concatenated_motionCorr_dfof'),:);
b2 = b_2.ctrlB2_meanWaveIntervals;
b2(isnan(b2))=[];

b_3 = neuronTable(strcmp(neuronTable.movie_name, '220215_n5_r1_concatenated_motionCorr_dfof'),:);
b3 = b_3.ctrlB2_meanWaveIntervals;
b3(isnan(b3))=[];

b_4 = neuronTable(strcmp(neuronTable.movie_name, '220224_n1_r1_concatenated_motionCorr_dfof'),:);
b4 = b_4.ctrlB2_meanWaveIntervals;
b4(isnan(b4))=[];

b_5 = neuronTable(strcmp(neuronTable.movie_name, '220225_n1_r2_concatenated_motionCorr_dfof'),:);
b5 = b_5.ctrlB2_meanWaveIntervals;
b5(isnan(b5))=[];

hex_iwi = neuronTable.hex_meanWaveIntervals;
hex_iwi(isnan(hex_iwi))=[];

h_1 = neuronTable(strcmp(neuronTable.movie_name, '220215_n2_r1_concatenated_motionCorr_dfof'),:);
h1 = h_1.hex_meanWaveIntervals;
h1(isnan(h1))=[];

h_2 = neuronTable(strcmp(neuronTable.movie_name, '220215_n4_r1_concatenated_motionCorr_dfof'),:);
h2 = h_2.hex_meanWaveIntervals;
h2(isnan(h2))=[];

h_3 = neuronTable(strcmp(neuronTable.movie_name, '220215_n5_r1_concatenated_motionCorr_dfof'),:);
h3 = h_3.hex_meanWaveIntervals;
h3(isnan(h3))=[];

h_4 = neuronTable(strcmp(neuronTable.movie_name, '220224_n1_r1_concatenated_motionCorr_dfof'),:);
h4 = h_4.hex_meanWaveIntervals;
h4(isnan(h4))=[];

h_5 = neuronTable(strcmp(neuronTable.movie_name, '220225_n1_r2_concatenated_motionCorr_dfof'),:);
h5 = h_5.hex_meanWaveIntervals;
h5(isnan(h5))=[];

mfa_iwi = neuronTable.mfa_meanWaveIntervals;
mfa_iwi(isnan(mfa_iwi))=[];

m_1 = neuronTable(strcmp(neuronTable.movie_name, '220215_n2_r1_concatenated_motionCorr_dfof'),:);
m1 = m_1.mfa_meanWaveIntervals;
m1(isnan(m1))=[];

m_2 = neuronTable(strcmp(neuronTable.movie_name, '220215_n4_r1_concatenated_motionCorr_dfof'),:);
m2 = m_2.mfa_meanWaveIntervals;
m2(isnan(m2))=[];

m_3 = neuronTable(strcmp(neuronTable.movie_name, '220215_n5_r1_concatenated_motionCorr_dfof'),:);
m3 = m_3.mfa_meanWaveIntervals;
m3(isnan(m3))=[];

m_4 = neuronTable(strcmp(neuronTable.movie_name, '220224_n1_r1_concatenated_motionCorr_dfof'),:);
m4 = m_4.mfa_meanWaveIntervals;
m4(isnan(m4))=[];

m_5 = neuronTable(strcmp(neuronTable.movie_name, '220225_n1_r2_concatenated_motionCorr_dfof'),:);
m5 = m_5.mfa_meanWaveIntervals;
m5(isnan(m5))=[];

%% 3) Testing distributions' normality and significance

[ctrl_H, ~,ctrl_swStat] = swtest(ctrl_iwi); %not normal H=1, reject the null
[ctrlB2_H, ~,ctrlB2_swStat] = swtest(ctrlB2_iwi); %not normal H=1, reject the null
[hex_H, ~,hex_swStat] = swtest(hex_iwi); %not normal H=1, reject the null
[mfa_H, ~,mfa_swStat] = swtest(mfa_iwi); %not normal H=1, reject the null

%two sample ks test compares the cdf of two data sets. 0 = cdfs are the
%same; 1 = cdfs are different.
[wt_b2_h, wt_b2_p,~] = kstest2(ctrl_iwi,ctrlB2_iwi); %wt vs ctrl
[b2C_b2H_h, b2C_b2H_p,~] = kstest2(ctrlB2_iwi,hex_iwi); %b2ctrl vs b2hex
[b2C_b2M_h, b2C_b2M_p,~] = kstest2(ctrlB2_iwi,mfa_iwi); %b2ctrl vs b2mfa
[b2M_b2H_h, b2M_b2H_p,~] = kstest2(hex_iwi,mfa_iwi); %b2mfa vs b2hex

b2C = kstest(ctrlB2_iwi);
% wilcoxon rank sum test to compare the medians of the wt and b2 ctrl
% distributions (for unparied data). H= 0 medians are the same, H = 1 medians are different.
% Gives two-sided p-value, doesn't test a specific null (i.e. median of wt
% is greater than median of b2) just looks at over diff in medians.
[WvB_p, WvB_h] = ranksum(ctrl_iwi,ctrlB2_iwi); 

% wilcoxon signed rank test to compare the medians of the b2 ctrl, hex and mfa
% distributions (for paried data). H= 0 medians are the same, H = 1 medians are different.
% Gives two-sided p-value, doesn't test a specific null (i.e. median of wt
% is greater than median of b2) just looks at over diff in medians.
data = [ctrl_iwi;ctrlB2_iwi;hex_iwi;mfa_iwi];

group1 = repmat("ctrl",length(ctrl_iwi),1);
group2 = repmat("ctrlB2",length(ctrlB2_iwi),1);
group3 = repmat("hex",length(hex_iwi),1);
group4 = repmat("mfa",length(mfa_iwi),1);
groups = [group1;group2;group3;group4];

[~,~,stats] = kruskalwallis(data, groups,"off");
[results,~,~,gnames] = multcompare(stats,"Alpha",0.01,"CriticalValueType","dunn-sidak");
tbl = array2table(results,"VariableNames", ...
    ["Group1","Group2","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group1") = gnames(tbl.("Group1"));
tbl.("Group2") = gnames(tbl.("Group2"))


%% 4) Plot data
% CDF of iwi 
path  = 'C:\Users\Christiane\Dropbox\Lab_Data\ipRGC_participation_in_waves\4. CompiledResults\Stage 1 waves\SpatioTemp\Summary figures\';

figure, cdfplot(ctrl_iwi)
hold on
cdfplot(ctrlB2_iwi)
cdfplot(hex_iwi)
cdfplot(mfa_iwi)
ylim([0 1])
xlim([0 700])
legend({'wt','b2','b2_hex','b2_mfa'})
savefig([path,'IWI_cdf_Summary.fig']);

% violin plots

figure
data = {ctrl_iwi, ctrlB2_iwi};
plotSpread(data, 'xNames', {'WT','B2_ctrl'}, 'distributionMarkers', {'.','.'},...
    'distributionColors', {'k','m'}, 'yLabel', 'mean iwi');
ylim([0 700])

figure
data1 = {hex_iwi, mfa_iwi};
plotSpread(data, 'xNames', {'B2_hex','B2_MFA'}, 'distributionMarkers', {'.','.'},...
    'distributionColors', {'c','y'}, 'yLabel', 'mean iwi');
ylim([0 700])

figure
violin({ctrl_iwi, ctrlB2_iwi, hex_iwi, mfa_iwi},'xlabel',{'WT','B2_ctrl','B2_hex','B2_MFA'},'facecolor',[0 0 0;1 0 1;0 1 1;1 1 0]);
ylim([0 700])

% plot and save iwi by condition

h1 = histogram(r1,edges, "Normalization","probability");
h1_values = h1.Values;

h2 = histogram(r2,edges, "Normalization","probability");
h2_values = h2.Values;

h3 = histogram(r3,edges, "Normalization","probability");
h3_values = h3.Values;

h4 = histogram(r4,edges, "Normalization","probability");
h4_values = h4.Values;

hb1 = histogram(b1,edges, "Normalization","probability");
hb1_values = hb1.Values;

hb2 = histogram(b2,edges, "Normalization","probability");
hb2_values = hb2.Values;

hb3 = histogram(b3,edges, "Normalization","probability");
hb3_values = hb3.Values;

hb4 = histogram(b4,edges, "Normalization","probability");
hb4_values = hb4.Values;

hb5 = histogram(b5,edges, "Normalization","probability");
hb5_values = hb5.Values;

figure
subplot(4,1,1)
histogram(ctrl_iwi,edges, "Normalization","probability")
title('wt ctrl iwi')
ylabel('% ROIs')
ylim([0 0.5])
xlabel('seconds')
subplot(4,1,2)
histogram(ctrlB2_iwi,edges, "Normalization","probability")
title('b2 ctrl iwi')
ylabel('% ROIs')
ylim([0 0.5])
xlabel('seconds')
subplot(4,1,3)
histogram(hex_iwi,edges, "Normalization","probability")
title('b2 hex iwi')
ylabel('% ROIs')
ylim([0 0.5])
xlabel('seconds')
subplot(4,1,4)
histogram(mfa_iwi,edges, "Normalization","probability")
title('b2 mfa iwi')
ylabel('% ROIs')
ylim([0 0.5])
xlabel('seconds')
% savefig([path,'wt_b2Drugs_hist_iwiSummary.fig']);
% hold off
% 
% pause()
% clear all, close all

%% extra plots

ctrl_meanIWI = nanmean(d2,1);
ctrl_d = mean(ctrl_meanIWI);
ctrl_err = std(nanstd(d2,0,1));

ctrlB2_meanIWI = nanmean(d3_c2,1);
ctrlB2_d = mean(ctrlB2_meanIWI);
ctrlB2_err = std(nanstd(d3_c2,0,1));

hex_meanIWI = nanmean(d3_h,1);
hex_d = mean(hex_meanIWI);
hex_err = std(nanstd(d3_h,0,1));

mfa_meanIWI = nanmean(d3_m,1);
mfa_d = mean(mfa_meanIWI);
mfa_err = std(nanstd(d3_m,0,1));

figure
subplot(1,2,1)
data = {ctrl_meanIWI, ctrlB2_meanIWI, hex_meanIWI, mfa_meanIWI};
plotSpread(data, 'xNames', {'WT','B2_ctrl','B2_hex','B2_MFA'}, 'distributionMarkers', {'o','o','o','o'},...
    'distributionColors', {'k','m','c','y'}, 'yLabel', 'mean iwi');
subplot(1,2,2)
y = [ctrl_d ctrlB2_d hex_d mfa_d];
x = categorical({'WT','B2_ctrl','B2_hex','B2_MFA'});
x = reordercats(x,{'WT','B2_ctrl','B2_hex','B2_MFA'});
bar(x,y);
ylim([0 400])
hold on
errorbar(x,y,[ctrl_d-ctrl_err ctrlB2_d-ctrlB2_err hex_d-hex_err mfa_d-mfa_err],[ctrl_d+ctrl_err ctrlB2_d+ctrlB2_err hex_d+hex_err mfa_d+mfa_err])

figure
subplot(1,2,1)
data = {ctrlB2_meanIWI, hex_meanIWI, mfa_meanIWI};
plotSpread(data, 'xNames', {'B2_ctrl','B2_hex','B2_MFA'}, 'distributionMarkers', {'o','o','o'},...
    'distributionColors', {'m','c','y'}, 'yLabel', '% Cells Participating');
ylim([0 400])
yticks(0:50:400)
subplot(1,2,2)
boxplot([ctrlB2_meanIWI', hex_meanIWI', mfa_meanIWI'],'Notch','on','Labels',{'B2_ctrl','B2_hex','B2_MFA'},'Whisker',1)
ylim([0 400])
yticks(0:50:400)
