%% Analyse av pasient nr 17
% Vi har to dataset tilgjengelig

load 20190117T145728_IQ_Sepsis-4min_traces;
Tmax1 = Tmax; Tmean1 = Tmean; Tmin1 = Tmin; Ts1 = Ts;
load 20190123T105641_IQ_Sepsis-4min_traces;
Tmax2 = Tmax; Tmean2 = Tmean; Tmin2 = Tmin; Ts2 = Ts;

%% Normal plot of both datasets

figure();clf;
subplot(2,2,1);plot(Tmean1.tED, Tmean1.velocity);grid;
title("Vmean dataset 20190117")
subplot(2,2,2);plot(Tmean1.tED, Tmean1.ART);grid;
title("Pmean dataset 20190117")
subplot(2,2,3);plot(Tmean2.tED, Tmean2.velocity);grid;
title("Vmean dataset 20190123")
subplot(2,2,4);plot(Tmean2.tED, Tmean2.ART);grid;
title("Pmean dataset 20190123")

%% Obtain frequency spectra (DFT)

figure();clf;

% ----------- 20190117 -----------

Fs1 = 1/(Tmean1.tED(2)-Tmean1.tED(1)); % Samplingsfrekvens
N1 = length(Tmean1.ART); %
NFFT1 = 2^nextpow2(N1); % Find number of bins needed for DFT (334 -> 512)
freqaxis1 = linspace(-Fs1/2, Fs1/2, NFFT1+1); % Get both side within nyquist freq: BW=Fs/2

% Velocity
V1_DC = mean(Tmean1.velocity); % Get DC component and subtract from signal
V1_Windowed = hamming(N_mean1).*(Tmean1.velocity-V1_DC); % Use hamming window
V1 = fftshift( abs( fft(  V1_Windowed, NFFT1+1 ) ) ) / N_mean1; % Get the mirrored and normalized DFT 
subplot(2,2,1);plot(freqaxis1, V1);grid;
title("Vmean Freq spectra - dataset 20190117")

% Pressure
P1_DC = mean(Tmean1.ART); % Get DC component and subtract from signal
P1_Windowed = hamming(N_mean1).*(Tmean1.ART-P1_DC); % Use hamming window
P1 = fftshift( abs( fft(  P1_Windowed, NFFT1+1 ) ) ) / N_mean1; % Get the mirrored and normalized DFT 
subplot(2,2,2);plot(freqaxis1, V1);grid;
title("Pmean Freq spectra - dataset 20190117")


% ----------- 20190123 -----------

Fs2 = 1/(Tmean2.tED(2)-Tmean2.tED(1)); % Samplingsfrekvens
N2 = length(Tmean2.ART);
NFFT2 = 2^nextpow2(N2); % Find number of bins needed for DFT (334 -> 512)
freqaxis2 = linspace(-Fs2/2, Fs2/2, NFFT2+1); % Get both side within nyquist freq: BW=Fs/2

% Velocity
V2_DC = mean(Tmean2.velocity); % Get DC component and subtract from signal
V2_Windowed = hamming(N_mean2).*(Tmean2.velocity-V2_DC); % Use hamming window
V2 = fftshift( abs( fft(  V2_Windowed, NFFT2+1 ) ) ) / N_mean2; % Get the mirrored and normalized DFT 
subplot(2,2,3);plot(freqaxis2, V2);grid;
title("Vmean Freq spectra - dataset 20190123")

% Pressure
P2_DC = mean(Tmean2.ART); % Get DC component and subtract from signal
P2_Windowed = hamming(N_mean2).*(Tmean2.ART-P2_DC); % Use hamming window
P2 = fftshift( abs( fft(  P2_Windowed, NFFT2+1 ) ) ) / N_mean2; % Get the mirrored and normalized DFT 
subplot(2,2,4);plot(freqaxis2, V2);grid;
title("Pmean Freq spectra - dataset 20190123")


%% Heartrate comparison

figure();clf;
subplot(2,1,1);plot(Tmean1.tED, Tmean1.HR);grid;
title("Heartrate dataset 20190117")
subplot(2,1,2);plot(Tmean2.tED, Tmean2.HR);grid;
title("Heartrate dataset 20190123")

%% Pulsatile Index (PI)

figure();clf;   
PI1 = (Tmax1.velocity-Tmin1.velocity)./Tmean1.velocity; % (max-min)/mean
subplot(2,2,1);plot(Tmax1.tED, Tmax1.velocity, Tmax1.velocity, Tmin1.velocity, Tmax1.tED, Tmean1.velocity);grid();
title("Velocities - dataset 20190117"); legend('Tmax','Tmin')
subplot(2,2,2);plot(Tmax1.tED, PI1);grid;
title("Pulsatile Index - dataset 20190117")

PI2 = (Tmax2.velocity-Tmin2.velocity)./Tmean2.velocity; % (max-min)/mean
subplot(2,2,3);plot(Tmax2.tED, Tmax2.velocity, Tmax2.velocity, Tmin2.velocity, Tmax2.tED, Tmean2.velocity);grid();
title("Velocities - dataset 20190123"); legend('Tmax','Tmin')
subplot(2,2,4);plot(Tmax2.tED, PI2);grid;
title("Pulsatile Index - dataset 20190123")