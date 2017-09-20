openfig('Cap_64QAM.fig');
hold on;
title('FER vs. Capacity of 64-QAM');
load('snr_FER_1_2_QAM64.mat');
plot(snr_FER,3,'r*');
load('snr_FER_2_3_A_QAM64.mat');
plot(snr_FER,(6*2/3),'r*');
load('snr_FER_3_4_A_QAM64.mat');
plot(snr_FER,4.5,'r*');
load('snr_FER_5_6_QAM64.mat');
plot(snr_FER,(6*5/6),'r*');

legend('64-QAM','R = 1/2', 'R = 2/3', 'R = 3/4', 'R = 5/6');
xlim([-5,30]);
ylim([0, 6.5]);