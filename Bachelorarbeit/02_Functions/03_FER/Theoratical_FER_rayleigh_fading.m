
%% Initialization
%addpath('Bachelorarbeit/cml/');
%CmlStartup;
clear all;
tic;

datapath1 = '/nas/ei/home/ga96jul/Bachelarbeit/Bachelorarbeit/03_Data/FER_AWGN_2310_01';
datapath2 = '/nas/ei/home/ga96jul/Bachelarbeit/Bachelorarbeit/03_Data/FER_AWGN_2310_001';
plotpath = '/nas/ei/home/ga96jul/Bachelarbeit/Bachelorarbeit/04_Plots/FER_AWGN_2310';

h = waitbar(0,'Calculating...');
n = 576;                                                                   % 576:96:2304
rate = (1/2);
ind = 0;
% R = k/n
snr_dB = 0:1:5;
cnt = 1;
Frame_errors = 0;
EsNo = 10.^(snr_dB/10);

%% Loop
%

for l = 1:length(snr_dB)
    pause(1);
frames = 10000000;  
%% Encoder
[H_rows, H_cols, P] = InitializeWiMaxLDPC(rate, n);                        % creating H-Matrix r x n

k = length( H_cols) - length(P);

for iterations = 1:frames 
data = (randi(1,1,k));

codeword = LdpcEncode(data, H_rows, P);                                    % codewords c

perm = randperm(length(codeword));

interleaver = intrlv(codeword,perm);
%% Mapping
QPSK = CreateConstellation( 'QPSK'); % QAM, PSK, FSK, etc. possible
symbols = [];
symbols = Modulate(interleaver,QPSK);


%%%%%%%%%%%%%%%%%%%%%%%% new
%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Split symbols into blocks

sym_temp = [];
variance = (2*EsNo(l));
amplitude = sqrt(EsNo(l));

sym_temp = symbols;
QPSK_A = QPSK;
block_sym = sym_temp;

%% channel
noise = (1/sqrt(variance))*(randn(size(block_sym(1,:)))+1i*randn(size(block_sym(1,:)))); 

r = amplitude * block_sym + noise;

receiv_sym = r;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
% %% Demapping
 sym_ll = Demod2D(receiv_sym, sqrt(EsNo(l))*QPSK_A , 1 );                                 % transforms received symbols into log-likelihoods
% 
 llr = Somap(sym_ll);                                                       % soft demapping
% 
deinterleaver = deintrlv(llr,perm);

% %% Decoder
 [output, errors] = MpDecode(-deinterleaver, H_rows, H_cols, 50, 0, 1, 1, data );     % decoder

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
%sem = semilogy(snr_dB,Frame_error_rate);
inter = linspace(0,5,5000000);
pFER = interp1(snr_dB,Frame_error_rate,inter);
semilogy(inter,pFER);
hold on;
xlim([0 5]);
grid on;
try
    
    snr_FER = find(pFER < 0.01);
    snr_FER = snr_FER(1)/1000000 + 0;
    plot(snr_FER, 0.01, 'r*');
catch 
    disp('No FER under 0.01');
end
try
    snr_FER_001 = find(pFER < 0.001);
    snr_FER_001 = snr_FER_001(1)/1000000 + 0;
    plot(snr_FER_001, 0.001, 'g*');
catch
    disp('No FER under 0.001');
end
print(plotpath);
save(datapath1,'Frame_error_rate');
%save('FER_1_2_QPSK_1810.fig','sem');
save(datapath2,'snr_FER');

toc;