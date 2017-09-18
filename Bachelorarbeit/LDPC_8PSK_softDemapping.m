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
snr_dB = 5;
cnt = 1;

%% Encoder
[H_rows, H_cols, P] = InitializeWiMaxLDPC(rate, n);                        % creating H-Matrix r x n

k = length( H_cols) - length(P);

data = round(rand(1,k));

codeword = LdpcEncode(data, H_rows, P);                                    % codewords c

%% Mapping
MPSK = CreateConstellation( 'PSK' , 8);                                    % QAM, PSK, FSK, etc. possible

symbols = Modulate(codeword,MPSK);

for u = 0:(length(codeword)/3)-1
    mapped_bit = codeword(1+(u*3):3+(u*3));
    mapped_bit_i = bi2de(mapped_bit);
    switch mapped_bit_i
        case 0
            level_bit(u+1) = 0;
        case 1
            level_bit(u+1) = 1;
        case 2
            level_bit(u+1) = 3;
        case 3
            level_bit(u+1) = 2;
        case 4
            level_bit(u+1) = 7;
        case 5
            level_bit(u+1) = 6;
        case 6
            level_bit(u+1) = 4;
        case 7
            level_bit(u+1) = 5;
    end
end

MPSK_symbol = exp(1i * (2*pi/8) * level_bit);

% channel

EsNo = 10^(snr_dB/10);

variance = 1/(2*EsNo);

noise = sqrt(variance)*(randn(size(symbols))+j*randn(size(symbols))); 

r = symbols + noise;  % symbols

%% Demapping
sym_ll = Demod2D(r, MPSK, EsNo);                                           % transforms received symbols into log-likelihoods

llr = Somap(sym_ll);                                                       % soft demapping

a = sqrt(2)/2;
x = [1 a+j*a j -a+j*a -1 -a-j*a -j a-j*a];
for s = 1:(n/3)
s1(s) = -EsNo*(abs(r(s)-x(1)).^2);
s2(s) = -EsNo*(abs(r(s)-x(2)).^2);
s3(s) = -EsNo*(abs(r(s)-x(3)).^2);
s4(s) = -EsNo*(abs(r(s)-x(4)).^2);
s5(s) = -EsNo*(abs(r(s)-x(5)).^2);
s6(s) = -EsNo*(abs(r(s)-x(6)).^2);
s7(s) = -EsNo*(abs(r(s)-x(7)).^2);
s8(s) = -EsNo*(abs(r(s)-x(8)).^2);

    
    L_y(1+(3*(s-1))) = log((exp(s5(s))+exp(s6(s))+exp(s7(s))+exp(s8(s)))/(exp(s1(s))+exp(s2(s))+exp(s3(s))+exp(s4(s))));
    L_y(2+(3*(s-1))) = log((exp(s3(s))+exp(s4(s))+exp(s7(s))+exp(s8(s)))/(exp(s1(s))+exp(s2(s))+exp(s5(s))+exp(s6(s))));
    L_y(3+(3*(s-1))) = log((exp(s3(s))+exp(s4(s))+exp(s5(s))+exp(s6(s)))/(exp(s1(s))+exp(s2(s))+exp(s7(s))+exp(s8(s))));
    

end

llr_new = (L_y);


%% Decoder
[output, errors] = MpDecode(-llr_new, H_rows, H_cols, 20, 0, 1, 1, data );     % decoder

