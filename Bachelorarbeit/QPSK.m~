%% Capacity calculation
% This script is intended to calculate to capacity of Gaussian Noise, QPSK
% and after user input any M-PSK that is desired. Everything is plotted for
% comparison.
%
% Input arguments:
% - number of samples N
% - SNR range in dB
% - Number of Points for MPSK
% - Arbitrary points for constellation 
%
% author : Kevin Li
% date : 06.09.2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%log changes

% 08.09.2017 
% -added nonstandard constellation points
% - note: constellation points must be normalized (affects amplitude ->
% snr). With points higher 1 we increase the SNR in calculation. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleanup and Initilization
clear all;

%% Inputs

prompt1 = 'Number Of Samples:  ';
N = input(prompt1);
prompt2 = 'SNR_minimum:  ';
snr_min = input(prompt2);
prompt3 = 'SNR_maximum:  ';
snr_max = input(prompt3);
prompt4 = 'How high of MPSK wants to be calculated?\n';
M = input(prompt4);
if ( floor(M) == M)
    modulation_scheme = M;
else
    disp('Error. Not valid input');
end
prompt5 = 'Input complex points as vector for nonstandard PSK ([x1, x2, x3...) with xi = x+jy. \n'; 
x_points = input(prompt5);

%% Gaussian Noise Capacity

snr_G = (snr_min:1:snr_max);                                               %SNR in dB
snr_power =10.^(snr_G/(10));                                              %SNR in W
h_Noise = log2(pi*exp(1));                                                 % Entropy of Noise with variance 1 (complex)
cap_Gaussian = [];

for power = 1:length(snr_G)

p = snr_power(power);

h_signal = log2(pi*exp(1+p));                                              % Entropy of received signal (complex)

cap_Gaussian(power) = log2(1+p);  %h_signal - h_Noise;                                  %h_signal - h_Noise;
end
figure;

plot(snr_G,cap_Gaussian);

hold on;
title('Capacity Calculation for Gaussian and PSK');
xlabel('SNR in dB');
ylabel('Capacity in bits/symbol');
ylim([0,8]);

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
Y = QPSK_symbols + Noise;                                                  % Y = X + N
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

%% M-PSK Capacity
% same as above with added modulation_scheme

Cap_MPSK = [];
for power = 1:length(snr)
    bits = round(modulation_scheme*(rand(N,1)));
    Amplitude_symbols = sqrt(snr_power(power));
    MPSK_symbols = Amplitude_symbols * exp(1i*(2*pi/modulation_scheme)*bits);
    Noise = (1/sqrt(2))*(randn(N,1)+j*randn(N,1));
    Y = MPSK_symbols + Noise;
    X = Amplitude_symbols * exp(1i*(2*pi/modulation_scheme)*(1:(modulation_scheme)));
    
    Y_t = transpose(Y);
    X_t = transpose(X);
    y = repmat(Y_t,modulation_scheme,1);
    x = repmat(X_t,1,N);
    s = -(abs(y-x)).^2;
    
    prob_Y = (1/(pi*modulation_scheme))*(sum(exp(s),1));
    h_Y = sum(-log2(prob_Y))/N;
    
    cap_MPSK(power) = (h_Y - h_Noise);
end
plot(snr,cap_MPSK);
desc = strcat(num2str(modulation_scheme), '-PSK');
legend('Gaussian','QPSK',desc);

%% Non standard PSK distribution
% 

% input x = [x1, x2, x3,..., xN]

%x_points = [1+j 1-j -1+j -1-j 3+j 3-j 3+3j 3-3j 1+3j 1-3j -1+3j -1-3j -3+3j -3-3j];

n = length(x_points);
h_N = log2(pi*exp(1));
scale = sum(abs(x_points).^2)/length(x_points);
x_points = x_points./sqrt(scale);

for l = 1:length(snr_power)
    level = round((n-1)*rand(1,N)+1);
    level_matrix = zeros(N,n);
for k = 1:N
    level_matrix(k,level(k)) = 1;
end

x_points = x_points./sqrt(sum(abs(x_points).^2));
x_points_t = transpose(x_points);
Amplitude = sqrt(snr_power(l));
x_out = Amplitude*level_matrix*x_points_t;

Noise = (1/sqrt(2))*(randn(N,1)+1i*randn(N,1));

Y = x_out + Noise;
Y_t = transpose(Y);
X = Amplitude * x_points_t;

y = repmat(Y_t,n,1);
x = repmat(X,1,N);

s = -((abs(y - x)).^2);
    
prob_Y = (1/(pi*n))*(sum(exp(s),1));
    
h_Y = sum(-log2(prob_Y))/N;

cap_PSK(l) = h_Y - h_N;
end

plot(snr_G,cap_PSK);
desc2 = strcat(num2str(n), '-Point-Modulation');
legend('Gaussian','QPSK',desc, desc2);
