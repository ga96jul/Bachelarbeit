%signal energy
%calculate the mean symbol energy
% Bit energy for each symbol
Eb   = 0;
EbNo = 10^(Eb/10);
% Noise variance   
No   = Eb/EbNo;  
symbol = 1:10000;
noise = sqrt(1/2)*(randn(1,length(symbol))+j*randn(1,length(symbol))); 
p = fft(noise)/length(symbol);

semilogy(symbol,(p));
xlabel('bins');
ylabel('power spectrum');
ylim([-0.000001, 1]);
