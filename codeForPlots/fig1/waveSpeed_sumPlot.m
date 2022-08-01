%This code:
% 1) takes the array of wave speeds copied from a separate excel
%file
% 2) runs stats on the data 
% 3) makes a summary figure.
%% concatenate speed data
l_waveVel = [174.1658075; 183.9436667; 192.190149; 162.2612779; 189.8368251; 237.3368694; 179.2431575; 165.6813312; 144.5218287; 182.0663497]; %Large Stage 1 waves
l2_waveVel = [77.90558644; 92.40253599; 98.1860635; 117.3380651;128.5237364; 87.97212726 ;92.40411967; 92.06444749; 103.246318; 103.8986282]; %Stage 2 waves
waveVel = horzcat(l_waveVel, l2_waveVel);

%% stats
[speed_h,speed_p,CI, Stats] = ttest2(l_waveVel,l2_waveVel);

%% summary plot
speedFig = figure;
makeFactorPlot(waveVel,speedFig)
xticklabels({'WT Control','B2 KO Control'})
% ylim([-1 1])
