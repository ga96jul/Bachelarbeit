%% Tx/Rx-Chain (Encoder/Decoder and Mapping/Demapping) for QPSK
%  Small script intended to get first look into the CM -Library and using
%  given functions. LDPC encoding was used here. 
%  Used functions:
%  - [H_rows, H_cols, P] = InitializeWiMaxLDPC( rate, number_of_symbols)
%  - LdpcEncode( data, H_rows, P)
%  - CreateConstellation( 'modulation_scheme', M);
%  - Modulate( codeword, Modulation_scheme)
%  - Demod2D( received_symbols, modulation_Scheme, EsNo)
%  - SoMap(loglikelihood)
%  - [output, errors] = MpDecode(-llr, H_rows, H_cols, iter, 0, 1, 1, data)
%
%  After restarting MATLAB CmlStartup has to be run once from the cml
%  folder!

%% Initialization
%addpath('Bachelorarbeit/cml/');
%CmlStartup;

n = 576;                                                                   % 576:96:2304
rate = 1/2;                                                                % R = k/n
snr_dB = 5.1;
cnt = 1;

%% Encoder
[H_rows, H_cols, P] = InitializeWiMaxLDPC(rate, n);                        % creating H-Matrix r x n

k = length( H_cols) - length(P);

data = round(rand(1,k));

codeword = LdpcEncode(data, H_rows, P);                                    % codewords c

%% Mapping
QPSK = CreateConstellation( 'QAM', 16);                                       % QAM, PSK, FSK, etc. possible

symbols = Modulate(codeword,QPSK);

% channel

EsNo = 10^(snr_dB/10);

variance = 1/(2*EsNo);

noise = sqrt(variance)*(randn(size(symbols))+j*randn(size(symbols))); 

r = symbols + noise;

%% Demapping
sym_ll = Demod2D(r, QPSK, EsNo);                                           % transforms received symbols into log-likelihoods

llr = Somap(sym_ll);                                                       % soft demapping
x = [1 j -j -1];
for s = 1:k/4
s1(s) = -EsNo*(abs(r(s)-x(1)).^2);
s2(s) = -EsNo*(abs(r(s)-x(2)).^2);
s3(s) = -EsNo*(abs(r(s)-x(3)).^2);
s4(s) = -EsNo*(abs(r(s)-x(4)).^2);


    L_y(1+(2*(s-1))) = log((exp(s1(s))+exp(s2(s)))/(exp(s3(s))+exp(s4(s))));
    L_y(2+(2*(s-1))) = log((exp(s3(s))+exp(s1(s)))/(exp(s4(s))+exp(s2(s))));

end

llr_new = (L_y);


%% Decoder
[output, errors] = MpDecode(llr_new, H_rows, H_cols, 50, 0, 1, 1, data );     % decoder
