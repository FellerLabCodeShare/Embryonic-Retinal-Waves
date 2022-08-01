%% Build the FOV table
% This code needs an excel spreadsheet with metadata in it.
% Right now, it will read 'cellParticipationMetaData_210312.xlsx'
% Every row in that excell spreadsheet is a field of view.
% the 'movie_nameSummary' column contains the name of the movie used to
% generate the set of variables being loaded here. The variables for each
% movie were also saved under the movie name.

%Steps:
% 1) load the data table
% 2) iterate through the rows of the table to load vairables from each
% movie to measure the mean event amplitude and area in control and DHBE
% conditions
% 3) Run stats and plot the data

%% 1) Load the table 
data_guide_name = 'cellParticipationMetaData_210312.xlsx';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'expDate','age','genotype', 'condition',...
    'movie_nameSummary','movie_nameVariables','movie_path','GridROIMask_name','GridROIMask_path','indicator','min_activity','expType','notes'},'string');

expTable = readtable(data_guide_name, opts');

[num_files dummy] = size(expTable);


%% 2) Here is the for loop that goes through every file
fileCounter = 1;

ctrl_table = expTable(strcmp(expTable.condition,'control'),:);
ctrl_meanEventArea = nan(length(ctrl_table.num_files), 1);
ctrl_meanEventAmp = nan(length(ctrl_table.num_files), 1);
for i = 1:length(ctrl_table.num_files)   % this will run through all files
    load(ctrl_table.movie_nameSummary(i))
    movie_name = ctrl_table.movie_nameSummary{i};
    movie_path = ctrl_table.movie_path{i};
    ctrl_meanEventArea(i,:) = mean(disc_percCellPart);
    ctrl_meanEventAmp(i,:) = mean(eventAmp);
end


dhbe_table = expTable(strcmp(expTable.condition,'dhbe'),:);
dhbe_meanEventArea = nan(length(dhbe_table.num_files), 1);
dhbe_meanEventAmp = nan(length(dhbe_table.num_files), 1);
for t = 1:length(dhbe_table.num_files)   % this will run through all files
    load(dhbe_table.movie_nameSummary(t))
    movie_name = dhbe_table.movie_nameSummary{t};
    movie_path = dhbe_table.movie_path{t};
    dhbe_meanEventArea(t,:) = mean(disc_percCellPart);
    dhbe_meanEventAmp(t,:) = mean(eventAmp);
end
%% 3) run stats and plot the results

%paried t-tests
[dhbeArea_h,dhbeArea_p] = ttest(ctrl_meanEventArea,dhbe_meanEventArea);
[dhbeAmp_h,dhbeAmp_p] = ttest(ctrl_meanEventAmp,dhbe_meanEventAmp);

%ctrl vs dhbe event area plot
eventArea_sum = [ctrl_meanEventArea dhbe_meanEventArea];
area = figure;
makeFactorPlot(eventArea_sum,area)
xticklabels({'Control','DHBE'})
ylim([0 100])

%ctrl vs dhbe evnt amp plot
eventAmp_sum = [ctrl_meanEventAmp dhbe_meanEventAmp];
amp = figure;
makeFactorPlot(eventAmp_sum,amp)
xticklabels({'Control','DHBE'})
ylim([0 6.2])

% %fold change plots
% ampChange = dhbe_meanEventAmp./ctrl_meanEventAmp;
% areaChange = dhbe_meanEventArea./ctrl_meanEventArea;
% 
% foldChangeData = {areaChange,ampChange};
% figure
% plotSpread(foldChangeData, 'xNames', {'Area_D','Amplitude_D'}, 'distributionMarkers', {'o'},...
%     'distributionColors', {'k'}, 'yLabel', 'Fold Change');
% ylim([0 1.5])

% 
% pause()
% close all

% save('waveCellPartTable.mat','waveCellPartTable');

