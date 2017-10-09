%% Frame Error Rate for QPSK with block fading
%  Calculate the frame error rate with LDPC Encode/Decoding done with the
%  CM-Library.


%% Initialization
%addpath('Bachelorarbeit/cml/');
%CmlStartup;
clear all;
tic;
h = waitbar(0,'Calculating...');
n = 2304;                                                                   % 576:96:2304
rate = (5/6);
ind = 0;
% R = k/n
snr_dB = 25:0.5:35;
cnt = 1;
Frame_errors = 0;
EsNo = 10.^(snr_dB/10);

%% Loop
%

for l = 1:length(snr_dB)
    pause(1);
frames = 10000;  
%% Encoder
[H_rows, H_cols, P] = InitializeWiMaxLDPC(rate, n);                        % creating H-Matrix r x n

k = length( H_cols) - length(P);

for iterations = 1:frames 
data = round(rand(1,k));

codeword = LdpcEncode(data, H_rows, P);                                    % codewords c

%% Mapping
QPSK = CreateConstellation( 'QPSK'); % QAM, PSK,save('snr_FER_1_2_QPSK_10_3.mat','snr_FER'); FSK, etc. possible
symbols = [];
symbols = Modulate(codeword,QPSK);


%%%%%%%%%%%%%%%%%%%%%%%% new
%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Split symbols into blocks

sym_temp = [];
variance = 1/(2*EsNo(l));
amplitude = sqrt(EsNo(l));
fading = (1/sqrt(2))*((randn(1))+1i*(randn(1)));


sym_temp = symbols;
QPSK_A = QPSK;
block_sym = [1 sym_temp];

%% channel
noise = (sqrt(variance))*(randn(size(block_sym(1,:)))+1i*randn(size(block_sym(1,:)))); 

r = block_sym*fading + noise;

est_fad = r(1);
est_sym = r/est_fad;
l_sym = length(est_sym);
real_sym = est_sym(2:l_sym);
receiv_sym = real_sym;

p = abs(est_fad);
E_fad = p^2*(p*exp(-((p^2)/2)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
fad_coeff = ones(1,l_sym-1)*est_fad;   
    
    
% %% Demapping
 sym_ll = Demod2D(receiv_sym, QPSK_A , EsNo(l)*E_fad);                                 % transforms received symbols into log-likelihoods
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
inter = linspace(25,35,10000);
pFER = interp1(snr_dB,Frame_error_rate,inter);
snr_FER = find(pFER < 0.01);
snr_FER = snr_FER(1)/1000 + 25;

save('FER_5_6_QPSK_fad_TFull.mat','Frame_error_rate');
save('FER_plot_5_6_QPSK_fad_TFull.fig','sem');
save('snr_FER_5_6_QPSK_fad_TFull.mat','snr_FER');

toc;