%% load data

load 20190123T105641_IQ_Sepsis-4min_traces;

N = length(Tmean.tED)
% dt (Ts) = 50ms 
fs = Ts.t(2)-Ts.t(1)


%% Remove offset (DC component)

N = length(Tmean.tED)
fshift  = (-N/2:1:(N/2-1))*fs;
figure();clf;

% ----------- PRESSURE -----------
% Remove DC component
P_DC = mean(Tmean.ART);
PmeanNew = Tmean.ART-P_DC;
% Do fft with hamming window and make it mirrored around 0 Hzs
P = fftshift( abs( fft(  hamming(N).*PmeanNew ) ) );
subplot(2,2,1);plot(Tmean.tED, PmeanNew);grid;
title("Pmean without DC component")
subplot(2,2,2);plot(fshift, P);grid;xlim([-0.6 0.6])
title("Absolute pressure freq spectra (shifted)")

% ----------- VELOCITY -----------
% Remove DC component
V_DC = mean(Tmean.velocity);
VmeanNew = Tmean.velocity - V_DC;
% Do fft with hamming window and make it mirrored around 0 Hzs
V = fftshift( abs( fft(  hamming(N).*VmeanNew ) ) );
subplot(2,2,3);plot(Tmean.tED, VmeanNew);grid;
title("Pmean without DC component")
subplot(2,2,4);plot(fshift, V);grid;xlim([-0.6 0.6])
title("Absolute velocity freq spectra (shifted)")


%% Pulsatility Index (PI)

figure();clf;
PI = (Tmax.ART-Tmin.ART)./Tmean.ART;
subplot(2,1,1);plot(Tmax.tED, Tmax.ART, Tmax.tED, Tmin.ART, Tmax.tED, Tmean.ART);grid();
subplot(2,1,2);plot(Tmax.tED, PI);grid;