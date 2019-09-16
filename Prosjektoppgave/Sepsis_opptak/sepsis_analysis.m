%% Dataset fra pasient med sepsis
% Det er uvist i hvilket stadie av sepsis pasienten befinner seg

load 20190117T145728_IQ_Sepsis-4min_traces;

Nmean = length(Tmean.tED)
Nraw = length(Ts.t)
% dt (Ts) = 50ms 
fs = Ts.t(2)-Ts.t(1)

%% Vanlig FFT
N = length(Tmean.tED)
f  = (0:(N-1))*fs;
figure();clf;

subplot(2,2,1);plot(Tmean.tED, Tmean.velocity);grid;
V = fft(Tmean.velocity);
subplot(2,2,2);plot(f, V);grid;

subplot(2,2,3);plot(Tmean.tED, Tmean.ART);grid;
P=fft(Tmean.ART); 
subplot(2,2,4);plot(f, P);grid;

%% Med FFTshift
N = length(Tmean.tED)
fshift  = (-N/2:1:(N/2-1))*fs;
figure();clf;

subplot(2,2,1);plot(Tmean.tED, Tmean.velocity);grid;
V = fft(Tmean.velocity);
V_shift = fftshift(V);
subplot(2,2,2);plot(fshift, V_shift);grid;

subplot(2,2,3);plot(Tmean.tED, Tmean.ART);grid;
P=fft(Tmean.ART); 
P_shift = fftshift(P);
subplot(2,2,4);plot(fshift, P_shift);grid;

%% Med Hamming vindu

N = length(Tmean.tED)
fshift  = (-N/2:1:(N/2-1))*fs;
figure();clf;

subplot(2,2,1);plot(Tmean.tED, Tmean.velocity);grid;
V = fft(hamming(N).*Tmean.velocity);
V_shift = fftshift(V);
subplot(2,2,2);plot(fshift, V_shift);grid;

subplot(2,2,3);plot(Tmean.tED, Tmean.ART);grid;
P=fft(hamming(N).*Tmean.ART); 
P_shift = fftshift(P);
subplot(2,2,4);plot(fshift, P_shift);grid;



%% Enkel MA filtrering

figure();clf;
filtered1 = movmean(Tmean.ART, 500);
subplot(2,1,1);plot(Tmean.tED, filtered1);grid;
filtered2 = movmean(Ts.ART, 5000);
subplot(2,1,2);plot(Ts.t, filtered2);grid;

%% MA frekvensanalyse for trend

figure();clf;
N = length(filtered1)
f  = (0:(N-1))*fs;
plot(f, fft(filtered1))
%plot(f, mag2db(abs(fft(filtered1))))

%% Med fftshift

figure();clf;

N = length(filtered1)
fshift  = (-N/2:1:(N/2-1))*fs;
x = fftshift(abs(fft(filtered1)));
plot(f, abs(fft(filtered1)))
%plot(fshift, mag2db(abs(fftshift(fft(filtered1)))))

%% Curvefitting all

figure();clf;
curvefit = fit( Tmean.tED, Tmean.ART, 'spline' )
plot(curvefit, Tmean.tED, Tmean.ART)

%% Curvefitting portion

x = Tmean.ART(1330:1458);
t = Tmean.tED(1330:1458);
figure();clf;
curvefit = fit( t, x, 'spline' )
plot(curvefit, t, x)