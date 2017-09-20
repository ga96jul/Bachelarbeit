%% Cleanup and Initilization
clear all;

%% Inputs

prompt1 = 'Number Of Samples:  ';
N = input(prompt1);
prompt2 = 'SNR_minimum:  ';
snr_min = input(prompt2);
prompt3 = 'SNR_maximum:  ';
snr_max = input(prompt3);

%% QPSK Capacity


h_Noise = log2(pi*exp(1));                                                              
snr = snr_min:1:snr_max;
snr_power = 10.^(snr/10);
cap_QPSK = [];


for power = 1:length(snr)
Amplitude_symbols = sqrt(snr_power(power));                                % Scaling of Amlitude through SNR  
bits = round(4*rand(N,1));                                                 % random bits created for phaseshift (between 0-4 int)
QPSK_symbols = Amplitude_symbols * exp(1i*(2*pi/4)*bits);                  % points shift by 90° to create QPSK
Noise = (1/sqrt(2))*(randn(N,1)+1i*randn(N,1));                            % add complex noise
Y   = QPSK_symbols + Noise;                                                  % Y = X + N
Y_t = transpose(Y);
X = Amplitude_symbols * exp(1i*(2*pi/4)*(1:4));                            % create (0,1),(0,-1),(-1,0),(1,0)
X_t = transpose(X);

%abs(y-xi)^2
x = repmat(X_t,1,N);                                                       % duplicate X symbols into a 4xN matrix
y = repmat(Y_t,4,1);                                                       % duplicate received symbols into a 4xN matrix
s = -(abs(y - x)).^2;                                                      % calculate difference for probability calculation

%h_signal
prob_Y = (1/(4*pi))*sum(exp(s),1);                                         % probability calculation
h_Y = (sum(-log2(prob_Y)))/N;                                              % entropy of received signal
cap_QPSK(power) = (h_Y - h_Noise);                                         % capacity
end
plot (snr,cap_QPSK);
hold on;
xlabel('SNR in dB');
ylabel('Bits per Symbol');
legend('QPSK');
grid on;