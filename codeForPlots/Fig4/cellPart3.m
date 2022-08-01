% Created by Alex and updated by Christiane
% This code loads the neuronTable matrix, which contains information (response amplitude during recording and whether or not it participated in a wave) on
% all ROIs across all FOVs, and measures the mean response amplitude during a wave for all
% ROIs across all FOVs. This code also measures the proportion of
% waves/events that each ROI participated in. ROIs are separated into two
% categories, RGCs and ipRGCs. 

% Steps:
% 1) load the neuronTable variable and separate out the control data.
% Further separate out the ipRGC containing data from the non-ipRGC
% containing data (based on tdTom fluorescence in the red channel during
% imaging)
% 2) itereate through FOVs and calculate the mean response amplitude for
% all RGCs and ipRGCs, and the proportion of waves RGCs and ipRGCs participated in.
% 3) Run stats
% 4) Plot the data



%% 1) Load the neuron table
load neuronTable.mat;
ctrl_neuronTable = neuronTable(strcmp(neuronTable.condition,'control'),:);%just the control FOVs
ipRGC_neuronTable = ctrl_neuronTable(strcmp(ctrl_neuronTable.iprgcID,"yes"),:); %just the ipRGCs
RGC_neuronTable = ctrl_neuronTable(strcmp(ctrl_neuronTable.iprgcID,"none"),:);  %just the RGCs

%% 2) control data
% Two useful vars here for the rest of the code
listFOVs = unique(ctrl_neuronTable.movie_name);  %List of unique FOVs
numFOVs = length(listFOVs);  %total num of FOVs

total_ipRGCs = length(ipRGC_neuronTable.neuronNum); %num of ipRGCs
total_RGCs = length(RGC_neuronTable.neuronNum); %num of RGCs

% Initializing some vars that we will fill in the for loop
ctrl_wavePartRGC = []; 
ctrl_wavePartIPRGC = [];

ctrl_meanAmpRGC = [];
ctrl_meanAmpIPRGC = [];


allIPRGCids = [];
allMaxCellAmps = [];


for i = 1:numFOVs
%     for i = 1:1   %This line for testing code

        %     generate tempTable of current FOV
        tempTable = ctrl_neuronTable(strcmp(ctrl_neuronTable.movie_name,listFOVs(i)),:);
        numWaves = tempTable.numWaves(i);  %total num of waves
        numCells = size(tempTable,1);
        %     Calc amp of cells when they participate in waves
        indWhenWavesHappen = tempTable.cellPart;
        indWhenWavesHappen(isnan(indWhenWavesHappen)) = 0;
        cellAmps = tempTable.maxCellAmp;
        cellAmps(isnan(cellAmps)) = 0;
        ampOfCellsDuringWaves = indWhenWavesHappen.*cellAmps;

        [ii,~,v] = find(ampOfCellsDuringWaves); %These two lines calc mean without counting zeros
        meanCellAmps = accumarray(ii,v,[],@mean);

        tempTable.meanCellAmps = meanCellAmps;

        allIPRGCids = [allIPRGCids;tempTable.iprgcID];
        allMaxCellAmps = [allMaxCellAmps;meanCellAmps];

        %     extract the ipRGC and rgc info to intialize variable with ipRGC and rgc avg part per FOV
        ctrlTable = tempTable(strcmp(tempTable.iprgcID,'none'),:);
        ipRGCTable = tempTable(strcmp(tempTable.iprgcID,'yes'),:);

        %     Calc percent waves each cell type participated in 
        wavePartRGC = mean(sum(ctrlTable.cellPart,2, 'omitnan')./numWaves); % average number of waves all RGCs within the FOV participated in
        wavePartIPRGC = mean(sum(ipRGCTable.cellPart,2, 'omitnan')./numWaves); % average number of waves all ipRGCs within the FOV participated in
        ctrl_wavePartRGC = [ctrl_wavePartRGC;wavePartRGC]; %saving RGC wave participation information for each FOV, of which there are 29.
        ctrl_wavePartIPRGC = [ctrl_wavePartIPRGC;wavePartIPRGC]; %saving ipRGC wave participation information for each FOV, of which there are 29.
        RGC_cellAmp = mean(ctrlTable.meanCellAmps); %average response amplitude of all RGCs within the FOV
        ipRGC_cellAmp = mean(ipRGCTable.meanCellAmps); %average response amplitude of all ipRGCs within the FOV
        ctrl_meanAmpRGC = [ctrl_meanAmpRGC;RGC_cellAmp]; %saving RGC mean response amplitude information for each FOV, of which there are 29.
        ctrl_meanAmpIPRGC = [ctrl_meanAmpIPRGC;ipRGC_cellAmp]; %saving ipRGC mean response amplitude information for each FOV, of which there are 29.
end


%% 3) stats for the cell part and event amp
% SW test of normality and anderson-darling test both tell us that these
% distributions are not normal
RGC_cellPart_normTest = swtest(ctrl_wavePartRGC);  %normal
ipRGC_cellPart_normTest = swtest(ctrl_wavePartIPRGC); %normal

RGC_eventAmp_normTest = swtest(ctrl_meanAmpRGC);  %normal
ipRGC_eventAmp_normTest = swtest(ctrl_meanAmpIPRGC); %normal


% unpaired ttest
[wavePart_h, wavePart_p] = ttest2(ctrl_wavePartRGC, ctrl_wavePartIPRGC);
[eventAmp_h, eventAmp_p] = ttest2(ctrl_meanAmpRGC, ctrl_meanAmpIPRGC);

%% 4) now plot a spreadplot and a boxplot of same data:
%cell participation mean
ipRGC_meanPart = mean(ctrl_wavePartIPRGC);
RGC_meanPart = mean(ctrl_wavePartRGC);
figure

subplot(1,2,1)
data = {ctrl_wavePartRGC, ctrl_wavePartIPRGC};
plotSpread(data, 'xNames', {'RGC','ipRGC'}, 'distributionMarkers', {'o','o'},...
    'distributionColors', {'k','r'}, 'yLabel', '% Cells Participating');
ylim([0 1])
yticks(0:0.5:1)
plot([ones(1,numFOVs);2*ones(1,numFOVs)], [ctrl_wavePartRGC';ctrl_wavePartIPRGC'],'k-')

subplot(1,2,2)
boxplot([ctrl_wavePartRGC, ctrl_wavePartIPRGC],'Notch','on','Labels',{'RGC','ipRGC'},'Whisker',1)
ylim([0 1])
yticks(0:0.5:1)
title('Compare Random Data from Different Distributions')





%now plot a spreadplot and a boxplot of same data:

figure
subplot(1,2,1)
data = {ctrl_meanAmpRGC, ctrl_meanAmpIPRGC};
plotSpread(data, 'xNames', {'RGC','ipRGC'}, 'distributionMarkers', {'.','.'},'distributionColors', {'k','r'},...
    'yLabel', 'mean cell amplitude');
ylim([0 3])
plot([ones(1,numFOVs);2*ones(1,numFOVs)], [ctrl_meanAmpRGC'; ctrl_meanAmpIPRGC'],'k-')

subplot(1,2,2)
boxplot([ctrl_meanAmpRGC, ctrl_meanAmpIPRGC],'Notch','on','Labels',{'RGC','ipRGC'},'Whisker',1)
ylim([0 3])
title('RGC vs ipRGC mean event amp per FOV')


iprgcIDs = strcmp(allIPRGCids,"yes");
rgcIDs = strcmp(allIPRGCids,"none");

allIPRGCmaxCellAmps = allMaxCellAmps(iprgcIDs);
allRGCmaxCellAmps = allMaxCellAmps(rgcIDs);

figure

subplot(1,3,1)
data = {allRGCmaxCellAmps, allIPRGCmaxCellAmps};
plotSpread(data, 'xNames', {'RGC','ipRGC'}, 'distributionMarkers', {'o','o'},...
    'distributionColors', {'k','r'}, 'yLabel', 'mean cell amplitude');
ylim([0 7])
subplot(1,3,2)
violin({allRGCmaxCellAmps,allIPRGCmaxCellAmps},'xlabel',{'RGC','ipRGC'},'facecolor',[1 1 0;0 1 0]);
ylim([0 7])

subplot(1,3,3)
% boxplot(allMaxCellAmps,allIPRGCids,'Notch','on','Labels',{'RGC','ipRGC'},'Whisker',1)
ylim([0 7])
