Fs= fs_mean1; % insert here your frequency sampling in Hz
L=length(Tmean1.ART); 
NFFT = 2^nextpow2(L);
V_DC = mean(Tmean1.ART);
Y  = fft( hamming(L).*(Tmean1.ART-V_DC),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

