openfig('Cap_64QAM.fig');
hold on;
title('FER (0.001) of QAM16 with Gaussian fading');
load('snr_FER_1_2_QAM64_001_2.mat');
plot(snr_FER,3,'b*');
load('snr_FER_2_3_QAM64_001_2.mat');
plot(snr_FER,(6*(2/3)),'b*');
load('snr_FER_3_4_QAM64_001_2.mat');
plot(snr_FER,6*(3/4),'b*');
load('snr_FER_5_6_QAM64_001_2.mat');
plot(snr_FER,6*(5/6),'b*');
grid on;

xlabel('SNR in dB');
ylabel('Bits per symbol');
legend('QAM-16', 'R = 1/2', 'R = 2/3', 'R = 3/4', 'R = 5/6');
xlim([0,30]);
ylim([0, 6.5]);