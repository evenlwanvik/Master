%% Extract parameters over time for two-element Windkessel model.

addpath('Master/Prosjektoppgave/Sepsis_opptak/patient17/data') 
load 20190117T145728_IQ_Sepsis-4min_traces;
date = '17.01.2019';
%load 20190123T105641_IQ_Sepsis-4min_traces;
%date = '23.01.2019';
name = append('Sepsis 4min_traces ', date)

%% Analyse measurements and connect end-points to avoid FFT leakage

delay=0.07;
tED=Tmean.tED-delay;
Ncycles=length(tED);

% Get data
n=140; % what cycle to extract
tIdx = find( Ts.t>tED(n) & Ts.t<tED(n+20)); % all samples within tED1 - tED2 window
t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx);
Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);



% Plot before
figure();clf;sgtitle(append('Removing discontinuities - measurements - ', date));
x_range = [t(1),t(end)];
subplot(2,2,1);plot(t,v);title('Velocity');xlim(x_range);grid();
subplot(2,2,2);plot(t,p);title('Pressure');xlim(x_range);grid();

% Plot hamming
v_windowed = hamming(N,'symmetric').*v;
p_windowed = hamming(N,'symmetric').*p;
subplot(2,2,3);plot(t,v_windowed);title('Velocity after hamming');xlim(x_range);grid();
subplot(2,2,4);plot(t,p_windowed);title('Pressure after hamming');xlim(x_range);grid();

%% Analyse the frequency spectra of the measurements before and after

f = (0:1/(N-1):1)*fs;
V = fft(v); V_hamming = fft(v_windowed); 
P = fft(p); P_hamming = fft(p_windowed); 

% plot all
figure();clf;sgtitle(append('V and P DFT - ', date));
subplot(2,2,1);plot(f,abs(V_hamming));title('Velocity');grid();xlim([0 500]);ylim([0, 0.1]);
subplot(2,2,2);plot(f,abs(P_hamming));title('Pressure');grid();xlim([0 500]);ylim([0, 10]);
subplot(2,2,3);plot(f,abs(V_after));title('Velocity Hamming');grid();xlim([0 500]);ylim([0, 0.1]);
subplot(2,2,4);plot(f,abs(P_after));title('Pressure Hamming');grid();xlim([0 500]);ylim([0, 10]);

%% Analyse the frequency spectra of the impedance before and after

w=2*pi*f;

% Get the index and 1. harmonic frequency
[Av, i_Fv] = max(abs(V(2:N-1))); i_Fv=i_Fv+1;

% Before
Z=P./V; 
R=Z(1);
C = abs( 1i*(1/Z(i_Fv)-1/R) / (2*pi*f(i_Fv)) )
%C = 1e-3;
Adm= 1/R - 1i*w*C;
Z_WK=1./Adm;

% After
Z_hamming=P_hamming./V_hamming;
R_hamming=Z_after(1);
C_hamming = abs( 1i*(1/Z_hamming(i_Fv)-1/R_hamming) / (2*pi*f(i_Fv)) );
Adm_WK_hamming= 1/R_hamming - 1i*w*C_hamming;
Z_WK_hamming=1./Adm_WK_hamming;

% plot
figure();clf;sgtitle(append('Z DFT - ', date));
subplot(2,1,1);plot(f,abs(Z),f,abs(Z_WK));title('Z before');grid();xlim([0,500]);ylim([0,2e4])
subplot(2,1,2);plot(f,abs(Z_hamming),f,abs(Z_WK_hamming));title('Z after');grid();ylim([0,2e4])


