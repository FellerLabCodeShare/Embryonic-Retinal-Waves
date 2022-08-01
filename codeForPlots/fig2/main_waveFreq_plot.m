 % created by Christiane with some code adapted from FSCH, previous graduate student in the lab (table2SubjectByConditionArray & makeFactorPlot)
 % This code loads the wave frequency meta data table and rearranges this
 % information to give you the frequency information for each FOV matching
 % each control FOV to its drug counterpart.

 % Steps:
 % 1) load the data table
 % 2) separate the data by experiment type/drug condition
 % 3) for each condition create an array with two columns, one for the
 % control data and the other for the drug data, making sure to keep the
 % information from the same FOV in the same row
 % 4) Plot data
 % 5) Run stats
%% 1) Load the table
data_guide_name = 'wavesFreqMetaData.xlsx';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
% opts = setvartype(opts,{'movie_name','wave_size', 'condition'},'string');

expTable = readtable(data_guide_name);

%% 2) creating sub_data tables based on experiment types
% ctrlExp = expTable(strcmp(expTable.expType,'control'),:);
% unique_ctrlCond = unique(ctrlExp.condition);
% schExp = expTable(strcmp(expTable.expType,'sch'),:);
% unique_schCond = unique(schExp.condition);
mfaExp = expTable(strcmp(expTable.expType,'mfa'),:);
unique_mfaCond = unique(mfaExp.condition);
dhbeExp = expTable(strcmp(expTable.expType,'dhbe'),:);
unique_dhbeCond = unique(dhbeExp.condition);
% lightExp = expTable(strcmp(expTable.expType,'light'),:);
% unique_lightCond = unique(lightExp.condition);
epiExp = expTable(strcmp(expTable.expType,'epi'),:);
unique_epiCond = unique(epiExp.condition);
hexExp = expTable(strcmp(expTable.expType,'hex'),:);
unique_hexCond = unique(hexExp.condition);

%% 3) new sub_data tables
% [ctrl_array, ctrl_subjects, ctrl_conditions] = table2SubjectByConditionArray(ctrlExp, unique_ctrlCond, ctrlExp.FOV_ID,  ctrlExp.condition, ctrlExp.waveFreq);
% [sch_array, sch_subjects, sch_conditions] = table2SubjectByConditionArray(schExp, unique_schCond, schExp.FOV_ID, schExp.condition, schExp.waveFreq);
[mfa_array, mfa_subjects, mfa_conditions] = table2SubjectByConditionArray(mfaExp, unique_mfaCond, mfaExp.FOV_ID, mfaExp.condition, mfaExp.waveFreq);
[dhbe_array, dhbe_subjects, dhbe_conditions] = table2SubjectByConditionArray(dhbeExp, unique_dhbeCond, dhbeExp.FOV_ID, dhbeExp.condition, dhbeExp.waveFreq);
% [light_array, light_subjects, light_conditions] = table2SubjectByConditionArray(lightExp, unique_lightCond, lightExp.FOV_ID, lightExp.condition, lightExp.waveFreq);
[epi_array, epi_subjects, epi_conditions] = table2SubjectByConditionArray(epiExp, unique_epiCond, epiExp.FOV_ID, epiExp.condition, epiExp.waveFreq);
[hex_array, hex_subjects, hex_conditions] = table2SubjectByConditionArray(hexExp, unique_hexCond, hexExp.FOV_ID, hexExp.condition, hexExp.waveFreq);

%% 4) making plots

% control plot
% ctrl = figure;
% makeFactorPlot(ctrl_array,ctrl)
% xticklabels({'Control'})
% ylim([0 2.5])

% sch plot
% sch = figure; 
% makeFactorPlot(sch_array,sch)
% xticklabels({'Control','SCH23390'})
% ylim([0 2.5])

%mfa plot
mfa = figure;
makeFactorPlot(mfa_array,mfa)
xticklabels({'Control','MFA'})
ylim([0 2])

%dhbe plot
% flip_dhbeArray = fliplr(dhbe_array);
dhbe = figure;
makeFactorPlot(dhbe_array,dhbe)
xticklabels({'Control','DHBE'})
ylim([0 2])

%epi plot
epi = figure;
makeFactorPlot(epi_array,epi)
xticklabels({'Control','Epi'})
ylim([0 2])

%hex plot
hex = figure;
makeFactorPlot(hex_array,hex)
xticklabels({'Control','Hex'})
ylim([0 1])

% %fold change plot
% fold_change = figure;
% % schFoldChange = sch_array(:,2)./sch_array(:,1);
% mfaFoldChange = mfa_array(:,2)./mfa_array(:,1);
% dhbeFoldChange = dhbe_array(:,2)./dhbe_array(:,1);
% epiFoldChange = epi_array(:,2)./epi_array(:,1);
% hexFoldChange = hex_array(:,2)./hex_array(:,1);
% 
% foldChangeData = {mfaFoldChange, dhbeFoldChange, epiFoldChange, hexFoldChange};
% plotSpread(foldChangeData, 'xNames', {'MFA','DHBE','Epi','Hex'}, 'distributionMarkers', {'o'},...
%     'distributionColors', {'k'}, 'yLabel', 'Fold Change');
% ylim([0 1.5])

%% 5) stats SW test and student T-test

% %sch data test for normality using the Shapiro-Wilk test. testing the distributions
% %individually
% sch_ctrl_normTest = swtest(sch_array(:,1)); %normal
% sch_sch_normTest = swtest(sch_array(:,2)); %normal
% % sch student paired t-test
% [sch_h,sch_p] = ttest(sch_array(:,1),sch_array(:,2));
% 
%mfa data test for normatlity
mfa_ctrl_normTest = swtest(mfa_array(:,1)); %normal
mfa_mfa_normTest = swtest(mfa_array(:,2)); %normal
% mfa student paired t-test
[mfa_h,mfa_p] = ttest(mfa_array(:,1),mfa_array(:,2));

%dhbe data test for normatlity

dhbe_ctrl_normTest = swtest(dhbe_array(:,1)); %normal
dhbe_dhbe_normTest = swtest(dhbe_array(:,2)); %normal
%dhbe student paired t-test (this is where doing a test to look at changes in shape
%of wave size distribution might great)
[dhbe_h,dhbe_p] = ttest(dhbe_array(:,1),dhbe_array(:,2));

%epi data test for normatlity
epi_ctrl_normTest = swtest(epi_array(:,1)); %normal
epi_epi_normTest = swtest(epi_array(:,2)); %not normal
%epi student paired t-test  
[epi_h,epi_p] = ttest(epi_array(:,1),epi_array(:,2));
epi_means = [mean(epi_array(:,1)), mean(epi_array(:,2))];
epi_stds = [std(epi_array(:,1)), std(epi_array(:,2))];

%hex data test for normatlity
hex_ctrl_normTest = swtest(hex_array(:,1)); %not normal
hex_hex_normTest = swtest(hex_array(:,2)); %normal
%hex student paired t-test 
[hex_h,hex_p] = ttest(hex_array(:,1),hex_array(:,2));

% %% fold change stats sw normality test and ANOVA
% 
% %normality test
% sch_foldChange_norm = swtest(schFoldChange); %normal
% mfa_foldChange_norm = swtest(mfaFoldChange); %normal
% dhbe_foldChange_norm = swtest(dhbeFoldChange); %normal
% epi_foldChange_norm = swtest(epiFoldChange); %not normal
% hex_foldChange_norm = swtest(hexFoldChange); %not normal
% 
% %anova
% data = padconcatenation(schFoldChange, mfaFoldChange,2);
% data_1 = padconcatenation(data,dhbeFoldChange,2);
% data_2 = padconcatenation(data_1, epiFoldChange,2);
% data_3 = padconcatenation(data_2, hexFoldChange,2);
% groups = {'SCH23390','MFA','dhbe','epi','hex'};
% [foldChange_p,~,stats] = anova1(data_3,groups);
% foldChange_stats = multcompare(stats);%default test for function is Tukey HSD test




