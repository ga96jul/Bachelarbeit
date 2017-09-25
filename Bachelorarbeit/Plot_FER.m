openfig('Cap_vs_FER_16QAM.fig');
hold on;
title('FER (0.001) of QAM-16 with Gaussian fading');
load('snr_FER_1_2_QAM16_001.mat');
plot(snr_FER,2,'b*');
load('snr_FER_2_3_QAM16_001.mat');
plot(snr_FER,(4*(2/3)),'b*');
load('snr_FER_3_4_QAM16_001.mat');
plot(snr_FER,3,'b*');
load('snr_FER_5_6_QAM16_001.mat');
plot(snr_FER,4*(5/6),'b*');
grid on;

xlabel('SNR in dB');
ylabel('Bits per symbol');
legend('QAM-16', 'R = 1/2', 'R = 2/3', 'R = 3/4', 'R = 5/6', 'In blue FER 0.001');
xlim([-5,25]);
ylim([0, 4.5]);