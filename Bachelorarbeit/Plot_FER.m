openfig('FER_vs_Capacity.fig');
hold on;
title('FER (0.001) of QPSK with Gaussian fading');
load('snr_FER_1_2_QPSK_001.mat');
plot(snr_FER,1,'b*');
load('snr_FER_2_3_QPSK_001.mat');
plot(snr_FER,(2*(2/3)),'b*');
load('snr_FER_3_4_QPSK_001.mat');
plot(snr_FER,1.5,'b*');
load('snr_FER_5_6_QPSK_001.mat');
plot(snr_FER,2*(5/6),'b*');
grid on;

xlabel('SNR in dB');
ylabel('Bits per symbol');
legend('QPSK', 'R = 1/2', 'R = 2/3', 'R = 3/4', 'R = 5/6', 'In blue FER 0.001');
xlim([-10,20]);
ylim([0, 2.5]);