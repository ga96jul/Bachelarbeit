%% soft decision decoding


x_points = [1+j 1-j -1+j -1-j ];
n = length(x_points);
N = 4;

snr_dB = 20;
snr_P = 10.^(snr_dB/10);
h_N = log2(pi*exp(1));

%for l = 1:length(snr_P)
    level = round((n-1)*rand(1,N)+1);
    level_matrix = zeros(N,n);
for k = 1:N
    level_matrix(k,level(k)) = 1;
end

x_points_t = transpose(x_points);
Amplitude = snr_P;
x_out = Amplitude*level_matrix*x_points_t;

Noise = (1/sqrt(2))*(randn(N,1)+j*randn(N,1));

Y = x_out + Noise;
Y_t = transpose(Y);

X = Amplitude * x_points_t;

%y = repmat(Y_t,n,1);
x = repmat(X,1,N);

for s = 1:N
s1(s) = -(abs(Y(s)-x(4)).^2);   %00
s2(s) = -(abs(Y(s)-x(3)).^2);   %10
s3(s) = -(abs(Y(s)-x(2)).^2);   %01
s4(s) = -(abs(Y(s)-x(1)).^2);   %11

    L_y(1+(2*(s-1))) = (s3(s)+log(1+exp(s1(s)-s3(s))))/(s4(s)+log(1+exp(s2(s)-s4(s))));
    L_y(2+(2*(s-1))) = (s3(s)+log(1+exp(s4(s)-s3(s))))/(s1(s)+log(1+exp(s2(s)-s1(s))));
end



%end