%% Frame Error Rate for QPSK with block fading
%  Calculate the frame error rate with LDPC Encode/Decoding done with the
%  CM-Library.


%% Initialization
%addpath('Bachelorarbeit/cml/');
%CmlStartup;
clear all;
tic;
h = waitbar(0,'Calculating...');noise = [];
n = 576;                                                                   % 576:96:2304
rate = (1/2);
ind = 0;
% R = k/n
snr_dB = 0:0.5:35;
cnt = 1;
Frame_errors = 0;
EsNo = 10.^(snr_dB/10);

T = n/16;
E_sum = 0;
noise = [];


%% Loop
%

for l = 1:length(snr_dB)
    pause(1);
frames = 10000;  
%% Encoder
[H_rows, H_cols, P] = InitializeWiMaxLDPC(rate, n);                        % creating H-Matrix r x n
%fad_coeff = ones(1,l_sym-1)*est_fad;  
k = length( H_cols) - length(P);

for iterations = 1:frames 
data = round(rand(1,k));

codeword = LdpcEncode(data, H_rows, P);                                    % codewords c

%% Mapping
QPSK = CreateConstellation( 'QPSK'); % QAM, PSK, FSK, etc. possible
symbols = [];
symbols = Modulate(codeword,QPSK);


%%%%%%%%%%%%%%%%%%%%%%%% new
%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Split symbols into blocks

sym_temp = [];
variance = 1/(2*EsNo(l));
amplitude = sqrt(EsNo(l));


for p = 1:1:(length(symbols)/T)

    
fading(p) = (1/sqrt(2))*((randn(1))+1i*(randn(1)));
sym_temp(p,:) = symbols((p-1)*T+1:p*T);
QPSK_A = QPSK;
block_sym(p,:) = [1 sym_temp(p,:)];

%% channel

noise(p,:) = (sqrt(variance))*(randn(size(block_sym(1,:)))+1i*randn(size(block_sym(1,:)))); 

r(p,:) = block_sym(p,:)*fading(p) + noise(p,:);

est_fad(p) = r(p,1);
est_sym(p,:) = r(p,:)/est_fad(p);
l_sym = length(est_sym(p,:));
real_sym(p,:) = est_sym(p,2:l_sym);

abs_f(p) = abs(est_fad(p));

end


receiv_sym = [];
receiv_sym = reshape(real_sym.',1,length(symbols));
sum_f = sum(abs_f);
for cnt = 1:1:length(symbols)/T
E_fad = abs_f(cnt)^2*(abs_f(cnt)*exp(-((abs_f(cnt)^2)/2)));
E_sum = E_sum + E_fad;
temp = E_sum;
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
% %% Demapping
 sym_ll = Demod2D(receiv_sym, QPSK_A , EsNo(l)*(E_sum/(length(symbols)/T)));    
 E_sum = 0;
 % transforms received symbols into log-likelihoods
% 
 llr = Somap(sym_ll);                                                       % soft demapping
% 
% %% Decoder
 [output, errors] = MpDecode(-llr, H_rows, H_cols, 50, 0, 1, 1, data );     % decoder

if (errors(50) > 0)
    Frame_errors = Frame_errors + 1;
end

if (Frame_errors == 100)
    break;
end

end
Frame_errors_SNR(l) = Frame_errors;
Frame_error_rate(l) = Frame_errors/iterations;
Frame_errors = 0;
waitbar(l/length(snr_dB));
end
close(h);
figure;
sem = semilogy(snr_dB,Frame_error_rate);
inter = linspace(20,35,15000);
pFER = interp1(snr_dB,Frame_error_rate,inter);
snr_FER = find(pFER < 0.01);
snr_FER = snr_FER(1)/1000 + 20;

%save('FER_1_2_QPSK_fad_T12.mat','Frame_error_rate');
%save('FER_plot_3_4_QPSK_fad_T12.fig','sem');
%save('snr_FER_1_2_QPSK_fad_T12.mat','snr_FER');

toc;