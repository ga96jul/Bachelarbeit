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
%% Non standard PSK distribution
% 

% input x = [x1, x2, x3,..., xN]

x_points = [-7-7j -7-5j -7-3j -7-1j -7+1j -7+3j -7+5j -7+7j -5-7j -5-5j -5-3j -5-1j -5+1j -5+3j -5+5j -5+7j -3-7j -3-5j -3-3j -3-1j -3+1j -3+3j -3+5j -3+7j -1-7j -1-5j -1-3j -1-1j -1+1j -1+3j -1+5j -1+7j 1-7j 1-5j 1-3j 1-1j 1+1j 1+3j 1+5j 1+7j 3-7j 3-5j 3-3j 3-1j 3+1j 3+3j 3+5j 3+7j 5-7j 5-5j 5-3j 5-1j 5+1j 5+3j 5+5j 5+7j 7-7j 7-5j 7-3j 7-1j 7+1j 7+3j 7+5j 7+7j];
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
figure;
plot(snr,cap_PSK);
desc = strcat('64-QAM');
hold on;
xlabel('SNR in dB');
ylabel('Bits per Symbol');
legend(desc);
grid on;