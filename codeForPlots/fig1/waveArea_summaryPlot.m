%% This code 
% 1) loads all wave area table (table with % of ROIs participating in each wave. since ROIs tiled the whole retina we used this percentage as a proxy for measuring wave size. 
% 2) Separate data by condition
% 3) Run the stats
% 4) plots a histogram with the mean value

%% 1) load wave area table
waveAreaVariable_path = 'C:\Users\Christiane\Dropbox\Lab_Data\ipRGC_participation_in_waves\4. CompiledResults\Stage 1 waves\SpatioTemp\Matlab variables for wave area';
cd(waveAreaVariable_path);
load('waveCellPartTable.mat');

%% 2) separate wt ctrl to b2 ctrl and all b2 conditions

%wave area
ctrl_table = waveCellPartTable(strcmp(waveCellPartTable.condition,'control'),:);
ctrl_waveArea = ctrl_table.ctrl_percCellPart;
ctrl_waveArea(isnan(ctrl_waveArea))=[];

r_1 = ctrl_table(strcmp(ctrl_table.movie_name, '211119_n3_r4_concatenated stack_motionCorr_dfof'),:);
r1 = r_1.ctrl_percCellPart;
r1(isnan(r1))=[];
r1 = r1';

r_2 = ctrl_table(strcmp(ctrl_table.movie_name, '211119_n4_r6_concatenated stack_motionCorr_dfof'),:);
r2 = r_2.ctrl_percCellPart;
r2(isnan(r2))=[];
r2 = r2';

r_3 = ctrl_table(strcmp(ctrl_table.movie_name, '211119_n6_r9_concatenated stack_motionCorr_dfof'),:);
r3 = r_3.ctrl_percCellPart;
r3(isnan(r3))=[];
r3 = r3';

r_4 = ctrl_table(strcmp(ctrl_table.movie_name, '211119_n6_r10_concatenated stack_motionCorr_dfof'),:);
r4 = r_4.ctrl_percCellPart;
r4(isnan(r4))=[];
r4 = r4';

b2_table = waveCellPartTable(strcmp(waveCellPartTable.condition,'b2ko_drug'),:);
ctrlB2_waveArea = b2_table.ctrlB2_percCellPart;
ctrlB2_waveArea(isnan(ctrlB2_waveArea))=[];

b_1 = b2_table(strcmp(b2_table.movie_name, '220215_n2_r1_concatenated_motionCorr_dfof'),:);
b1 = b_1.ctrlB2_percCellPart;
b1(isnan(b1))=[];
b1 = b1';

b_2 = b2_table(strcmp(b2_table.movie_name, '220215_n4_r1_concatenated_motionCorr_dfof'),:);
b2 = b_2.ctrlB2_percCellPart;
b2(isnan(b2))=[];
b2 = b2';

b_3 = b2_table(strcmp(b2_table.movie_name, '220215_n5_r1_concatenated_motionCorr_dfof'),:);
b3 = b_3.ctrlB2_percCellPart;
b3(isnan(b3))=[];
b3 = b3';

b_4 = b2_table(strcmp(b2_table.movie_name, '220224_n1_r1_concatenated_motionCorr_dfof'),:);
b4 = b_4.ctrlB2_percCellPart;
b4(isnan(b4))=[];
b4 = b4';

b_5 = b2_table(strcmp(b2_table.movie_name, '220225_n1_r2_concatenated_motionCorr_dfof'),:);
b5 = b_5.ctrlB2_percCellPart;
b5(isnan(b5))=[];
b5 = b5';

hex_waveArea = b2_table.hex_percCellPart;
hex_waveArea(isnan(hex_waveArea))=[];

h_1 = b2_table(strcmp(b2_table.movie_name, '220215_n2_r1_concatenated_motionCorr_dfof'),:);
h1 = h_1.hex_percCellPart;
h1(isnan(h1))=[];
h1 = h1';

h_2 = b2_table(strcmp(b2_table.movie_name, '220215_n4_r1_concatenated_motionCorr_dfof'),:);
h2 = h_2.hex_percCellPart;
h2(isnan(h2))=[];
h2 = h2';

h_3 = b2_table(strcmp(b2_table.movie_name, '220215_n5_r1_concatenated_motionCorr_dfof'),:);
h3 = h_3.hex_percCellPart;
h3(isnan(h3))=[];
h3 = h3';

h_4 = b2_table(strcmp(b2_table.movie_name, '220224_n1_r1_concatenated_motionCorr_dfof'),:);
h4 = h_4.hex_percCellPart;
h4(isnan(h4))=[];
h4 = h4';

h_5 = b2_table(strcmp(b2_table.movie_name, '220225_n1_r2_concatenated_motionCorr_dfof'),:);
h5 = h_5.hex_percCellPart;
h5(isnan(h5))=[];
h5 = h5';

mfa_waveArea = b2_table.mfa_percCellPart;
mfa_waveArea(isnan(mfa_waveArea))=[];

m_1 = b2_table(strcmp(b2_table.movie_name, '220215_n2_r1_concatenated_motionCorr_dfof'),:);
m1 = m_1.mfa_percCellPart;
m1(isnan(m1))=[];
m1 = m1';

m_2 = b2_table(strcmp(b2_table.movie_name, '220215_n4_r1_concatenated_motionCorr_dfof'),:);
m2 = m_2.mfa_percCellPart;
m2(isnan(m2))=[];
m2 = m2';

m_3 = b2_table(strcmp(b2_table.movie_name, '220215_n5_r1_concatenated_motionCorr_dfof'),:);
m3 = m_3.mfa_percCellPart;
m3(isnan(m3))=[];
m3 = m3';

m_4 = b2_table(strcmp(b2_table.movie_name, '220224_n1_r1_concatenated_motionCorr_dfof'),:);
m4 = m_4.mfa_percCellPart;
m4(isnan(m4))=[];
m4 = m4';

m_5 = b2_table(strcmp(b2_table.movie_name, '220225_n1_r2_concatenated_motionCorr_dfof'),:);
m5 = m_5.mfa_percCellPart;
m5(isnan(m5))=[];
m5 = m5';

% %activity index code that was not used
% ctrl_actIdx = ctrl_table.activity_idx;
% ctrlB2_actIdx = b2_table.ctrlB2Activity_idx;
% hex_actIdx = b2_table.hexActivity_idx;
% mfa_actIdx = b2_table.mfaActivity_idx;
% 
% c_data = padconcatenation(ctrl_actIdx,ctrlB2_actIdx,2);
% b2_data = horzcat(ctrlB2_actIdx,hex_actIdx,mfa_actIdx);
% control = figure;
% makeFactorPlot(c_data,control)
% xticklabels({'WT Control','B2 KO Control'})
% ylim([0 2])
% 
% b2_ko = figure;
% makeFactorPlot(b2_data,b2_ko)
% xticklabels({'B2 KO Control','B2 KO Hex','B2 KO MFA'})
% ylim([0 2])
%% 3) Stats

D = [ctrl_waveArea,ctrlB2_waveArea,hex_waveArea,mfa_waveArea];
group1 = repmat("ctrl",length(ctrl_waveArea),1);
group2 = repmat("ctrlB2",length(ctrlB2_waveArea),1);
group3 = repmat("hex",length(hex_waveArea),1);
group4 = repmat("mfa",length(mfa_waveArea),1);
G = [group1;group2;group3;group4];
G = G';
[~,~,stats] = anova1(D,G,"off");
[results,~,~,gnames] = multcompare(stats,'Alpha', 0.01,'CriticalValueType','tukey-kramer');
tbl = array2table(results,"VariableNames", ...
    ["Group1","Group2","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group1") = gnames(tbl.("Group1"));
tbl.("Group2") = gnames(tbl.("Group2"))


%% 4) Plot data 
% CDF of wave area
path  = 'C:\Users\Christiane\Dropbox\Lab_Data\ipRGC_participation_in_waves\4. CompiledResults\Stage 1 waves\SpatioTemp\Summary figures\';

figure, cdfplot(ctrl_waveArea)
hold on
cdfplot(ctrlB2_waveArea)
cdfplot(hex_waveArea)
cdfplot(mfa_waveArea)
ylim([0 1])
xlim([0 100])
legend({'wt','b2','b2_hex','b2_mfa'})
savefig([path,'waveArea_cdf_Summary.fig']);

% violin plot

figure
subplot(1,2,1)
data = {ctrl_waveArea, ctrlB2_waveArea, hex_waveArea, mfa_waveArea};
plotSpread(data, 'xNames', {'WT','B2_ctrl','B2_hex','B2_MFA'}, 'distributionMarkers', {'.','.','.','.'},...
    'distributionColors', {'k','m','c','y'}, 'yLabel', '% Retina');
ylim([0 100])
subplot(1,2,2)
violin({ctrl_waveArea, ctrlB2_waveArea, hex_waveArea, mfa_waveArea},'xlabel',{'WT','B2_ctrl','B2_hex','B2_MFA'},'facecolor',[0 0 0;1 0 1;0 1 1;1 1 0]);
ylim([0 100])

% plot wave area
path  = 'C:\Users\Christiane\Dropbox\Lab_Data\ipRGC_participation_in_waves\4. CompiledResults\Stage 1 waves\SpatioTemp\Summary figures\';

figure, histogram(ctrl_waveArea,0:10:100, "Normalization","probability")
title('wt wave area')
ylabel('% Waves')
xlabel('% Retina')
savefig([path,'wt_waveAreaSummary.fig']);

figure
subplot(4,1,1)
histogram(ctrl_waveArea,0:10:100, "Normalization","probability")
title('wt wave area')
ylabel('% Waves')
xlabel('% Retina')
ylim([0 0.5])
subplot(4,1,2)
histogram(ctrlB2_waveArea,0:10:100, "Normalization","probability")
title('b2 ctrl wave area')
ylabel('% Waves')
xlabel('% Retina')
ylim([0 0.5])
subplot(4,1,3)
histogram(hex_waveArea,0:10:100, "Normalization","probability")
title('b2 hex wave area')
ylabel('% Waves')
xlabel('% Retina')
ylim([0 0.5])
subplot(4,1,4)
histogram(mfa_waveArea,0:10:100, "Normalization","probability")
title('b2 mfa wave area')
ylabel('% Waves')
xlabel('% Retina')
ylim([0 1])
savefig([path,'b2_waveAreaDrugSummary.fig']);

% pause()
% clear all, close all