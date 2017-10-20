%% Non standard PSK distribution
% 

% input x = [x1, x2, x3,..., xN]

x_points = [j -1 -j 1 2 -2];
n = length(x_points);
N = 100000;

snr_dB = 1:30;
snr_P = 10.^(snr_dB/10);
h_N = log2(pi*exp(1));

for l = 1:length(snr_P)
    level = round((n-1)*rand(1,N)+1);
    level_matrix = zeros(N,n);
for k = 1:N
    level_matrix(k,level(k)) = 1;
end

x_points_t = transpose(x_points);
Amplitude = snr_P(l);
x_out = Amplitude*level_matrix*x_points_t;

Noise = (1/sqrt(2))*(randn(N,1)+j*randn(N,1));

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

plot(snr_dB,cap_PSK);