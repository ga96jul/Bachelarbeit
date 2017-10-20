load('FER_AWGN_1910');
SNR_dB = 1:50;
FER_ray = [];
for l = 1:length(SNR_dB)
    h = (1/sqrt(2))*((randn(1,100000))+1i*(randn(1,100000)));
    snr = 10.^(SNR_dB/10);
    power = abs(h).^2;
    snr_ray = power*SNR_dB(l);
    for p = 1:length(snr_ray)
        power = abs(h(p))^2;
        snr_rayleigh = power*snr(l);
        snr_real = 10*log10(snr_rayleigh);
        slot_FER = snr_real*10000;
        rounded = fix(slot_FER);
        try
            FER_ray(p) = pFER(rounded);
        catch
            FER_ray(p) = 0;
        end
    end
    FER_rayleigh(l) = sum(FER_ray)/100000;
    if (FER_rayleigh(l) > 1)
        FER_rayleigh(l) = 1;
    end
    
    
end

figure;
semilogy(SNR_dB,FER_rayleigh);