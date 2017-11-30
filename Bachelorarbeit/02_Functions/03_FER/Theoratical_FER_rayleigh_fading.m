
%% Initialization
%addpath('Bachelorarbeit/cml/');
%CmlStartup;
clear all;
tic;


h = waitbar(0,'Calculating...');
n = 576;                                                                   % 576:96:2304
rate = (1/2);
ind = 0;
% R = k/n
snr_dB = 0:1:6;
cnt = 1;
Frame_errors = 0;
EsNo = 10.^(snr_dB/10);

%% Loop
%

for l = 1:length(snr_dB)
    pause(1);
frames = 100000000;  
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
noise = (1/sqrt(2))*(randn(size(block_sym(1,:)))+1i*randn(size(block_sym(1,:)))); 

r = amplitude * block_sym + noise;

receiv_sym = r;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    


x = amplitude*[1 j -j -1];
temp_snr = EsNo(l);
for s = 1:n/2
s1(s) = -(1/EsNo(l))*(abs(receiv_sym(s)-x(1)).^2);    %1
s2(s) = -(1/EsNo(l))*(abs(receiv_sym(s)-x(2)).^2);    
s3(s) = -(1/EsNo(l))*(abs(receiv_sym(s)-x(3)).^2);
s4(s) = -(1/EsNo(l))*(abs(receiv_sym(s)-x(4)).^2);

if(s2(s)<s1(s))
    if(s4(s)<s3(s))

        L_y(1+(2*(s-1))) = (s1(s)+log(1+exp(s2(s)-s1(s))))-(s3(s)+log(1+exp(s4(s)-s3(s))));
    else
        L_y(1+(2*(s-1))) = (s1(s)+log(1+exp(s2(s)-s1(s))))-(s4(s)+log(1+exp(s3(s)-s4(s))));
    end
else
    if(s4(s)<s3(s))

        L_y(1+(2*(s-1))) = (s2(s)+log(1+exp(s1(s)-s2(s))))-(s3(s)+log(1+exp(s4(s)-s3(s))));
    else
        L_y(1+(2*(s-1))) = (s2(s)+log(1+exp(s1(s)-s2(s))))-(s4(s)+log(1+exp(s3(s)-s4(s))));
    end
end
  
if(s1(s)<s3(s))
    if(s2(s)<s4(s))
        L_y(2+(2*(s-1))) = (s3(s)+log((1+exp(s1(s)-s3(s)))))-(s4(s)+log(1+exp(s2(s)-s4(s))));
    else
        L_y(2+(2*(s-1))) = (s3(s)+log((1+exp(s1(s)-s3(s)))))-(s2(s)+log(1+exp(s4(s)-s2(s))));
    end
else
    if(s2(s)<s4(s))
        L_y(2+(2*(s-1))) = (s1(s)+log((1+exp(s3(s)-s1(s)))))-(s4(s)+log(1+exp(s2(s)-s4(s))));
    else
        L_y(2+(2*(s-1))) = (s1(s)+log((1+exp(s3(s)-s1(s)))))-(s2(s)+log(1+exp(s4(s)-s2(s))));
    end
end

%      L_y(1+(2*(s-1))) = (s1(s)+log(1+exp(s2(s)-s1(s))))/(s3(s)+log(1+exp(s4(s)-s3(s))));
%      L_y(2+(2*(s-1))) = (s3(s)+log(1+exp(s1(s)-s3(s))))/(s4(s)+log(1+exp(s2(s)-s4(s))));


% L_y(1+(4*(s-1))) = s1(s);
% L_y(2+(4*(s-1))) = s2(s);
% L_y(3+(4*(s-1))) = s3(s);
% L_y(4+(4*(s-1))) = s4(s);
end

llr_new = -L_y;



deinterleaver = deintrlv(llr_new,perm);

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
temp = sprintf('Frame_errors_%d',l);
temp2 = sprintf('FER_%d',l);
save(temp,'Frame_errors_SNR');
save(temp2,'Frame_error_rate');
end
close(h);
figure;
%sem = semilogy(snr_dB,Frame_error_rate);
inter = linspace(0,5,50000000);
pFER = interp1(snr_dB,Frame_error_rate,inter);
semilogy(inter,pFER);
hold on;
xlim([0 5]);
grid on;
try
    
    snr_FER = find(pFER < 0.01);
    snr_FER = snr_FER(1)/10000000 + 0;
    plot(snr_FER, 0.01, 'r*');
catch 
    disp('No FER under 0.01');
end
try
    snr_FER_001 = find(pFER < 0.001);
    snr_FER_001 = snr_FER_001(1)/10000000 + 0;
    plot(snr_FER_001, 0.001, 'g*');
catch
    disp('No FER under 0.001');
end
print(plotpath);
save(datapath1,'Frame_error_rate');
%save('FER_1_2_QPSK_1810.fig','sem');
save(datapath2,'snr_FER');

toc;