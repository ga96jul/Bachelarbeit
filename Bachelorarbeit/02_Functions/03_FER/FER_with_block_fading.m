%% Frame Error Rate for QPSK with block fading
%  Calculate the frame error rate with LDPC Encode/Decoding done with the
%  CM-Library.


%% Initialization
%addpath('Bachelorarbeit/cml/');
%CmlStartup;
clear all;
%estimate unknown fading
tic;
h = waitbar(0,'Calculating...');
n = 2304;                                                                   % 576:96:2304
rate = (1/2);
ind = 0;
% R = k/n
snr_dB = 10:0.5:15;
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
QPSK = CreateConstellation( 'QPSK');                                       % QAM, PSK, FSK, etc. possible
symbols = [];
symbols = Modulate(codeword,QPSK);

%% Split symbols into blocks

T = 4;
sym_temp = [];
variance = 1/(2*EsNo(l));
j = 0;
for j = 1:1:(length(symbols)/T)
    
    fading = (1/sqrt(2))*(randn(1)+1i*randn(1));
    sym_temp(j,:) = symbols((j-1)*T+1:j*T);
    block_sym(j,:) = [1 sym_temp(j,:)];

%% channel
    noise = sqrt(variance)*(randn(size(block_sym(1,:)))+1i*randn(size(block_sym(1,:)))); 

    r(j,:) = block_sym(j,:)*fading + noise;


    est_fad = r(j,1);
    est_sym = r(j,:)/est_fad;
    l_sym = size(est_sym);
    
    real_sym(j,:) = est_sym(2:l_sym(2));
    end
    receiv_sym = [];
    receiv_sym = reshape(real_sym.',1,length(symbols));

% %% Demapping
 sym_ll = Demod2D(receiv_sym, QPSK , EsNo*est_fad^2);       %adjust to fading                          % transforms received symbols into log-likelihoods
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
inter = linspace(20,25,5000);
pFER = interp1(snr_dB,Frame_error_rate,inter);
snr_FER = find(pFER < 0.01);
snr_FER = snr_FER(1)/1000 + 20;

save('FER_1_2_QPSK_Fad.mat','Frame_error_rate');
save('FER_plot_1_2_QPSK_Fad.fig','sem');
save('snr_FER_1_2_QPSK_Fad.mat','snr_FER');

toc;