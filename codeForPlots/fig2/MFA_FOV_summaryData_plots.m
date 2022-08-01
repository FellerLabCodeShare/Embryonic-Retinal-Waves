%% Build the FOV table

% This code needs an excel spreadsheet with metadata in it.
% Right now, it will read 'cellParticipationMetaDataMFA.xlsx'
% Every row in that excell spreadsheet is a field of view.
% the 'movie_nameSummary' column contains the name of the movie used to
% generate the set of variables being loaded here. The variables for each
% movie were also saved under the movie name.

%Steps:
% 1) load the data table
% 2) iterate through the rows of the table to load vairables from each
% movie to measure the mean event amplitude and area in control and MFA
% conditions
% 3) Run stats and plot the data

%% 1) Load the table 
data_guide_name = 'cellParticipationMetaData_MFA.xlsx';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'expDate','age','genotype', 'condition',...
    'movie_nameSummary','movie_nameVariables','movie_path','ROIMask_name','ROIMask_path','indicator','min_activity','expType','notes'},'string');

expTable = readtable(data_guide_name, opts');

[num_files dummy] = size(expTable);

%% 2) Here is the for loop that goes through every file
fileCounter = 1;

ctrl_table = expTable(strcmp(expTable.condition,'control'),:);
ctrl_meanEventArea = nan(length(ctrl_table.num_files), 1);
ctrl_meanEventAmp = nan(length(ctrl_table.num_files), 1);
for i = 1:length(ctrl_table.num_files)   % this will run through all files and load the set of matlab variables saved using the mive name for the file currently being observed
    load(ctrl_table.movie_nameSummary(i))
    movie_name = ctrl_table.movie_nameSummary{i};
    movie_path = ctrl_table.movie_path{i};
    ctrl_meanEventArea(i,:) = mean(disc_percCellPart);
    ctrl_meanEventAmp(i,:) = mean(eventAmp);
end


mfa_table = expTable(strcmp(expTable.condition,'mfa'),:);
mfa_meanEventArea = nan(length(mfa_table.num_files), 1);
mfa_meanEventAmp = nan(length(mfa_table.num_files), 1);
for t = 1:length(mfa_table.num_files)   % this will run through all files
    load(mfa_table.movie_nameSummary(t))
    movie_name = mfa_table.movie_nameSummary{t};
    movie_path = mfa_table.movie_path{t};
    mfa_meanEventArea(t,:) = mean(disc_percCellPart);
    mfa_meanEventAmp(t,:) = mean(eventAmp);
end
%% 3) run stats and plot the results

%paried t-tests
[mfaArea_h,mfaArea_p] = ttest(ctrl_meanEventArea,mfa_meanEventArea);
[mfaAmp_h,mfaAmp_p] = ttest(ctrl_meanEventAmp,mfa_meanEventAmp);


%ctrl vs dhbe event area plot
eventArea_sum = [ctrl_meanEventArea mfa_meanEventArea];
area = figure;
makeFactorPlot(eventArea_sum,area)
xticklabels({'Control','MFA'})
ylim([0 100])


%ctrl vs dhbe evnt amp plot
eventAmp_sum = [ctrl_meanEventAmp mfa_meanEventAmp];
amp = figure;
makeFactorPlot(eventAmp_sum,amp)
xticklabels({'Control','MFA'})
ylim([0 3.05])


% 
% pause()
% close all

% save('waveCellPartTable.mat','waveCellPartTable');

