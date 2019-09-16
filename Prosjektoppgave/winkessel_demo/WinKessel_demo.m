%% Winkessel demo
% Endre motstand R og kondensator C, og observer endring i flow kurveform
% Regn ut pulsativ index  (PI) for p og q
% Hvordan endrer PI seg med varierende R og C?

% 2019.08.19 . Hans Torp

load pressureMeasurement;% trykk p og tidsvektor t
figure(1);clf;
subplot(2,1,1);plot(t,p);grid;
title("Pressure")
%%
N=length(t);
dt=(t(end)-t(1))/(N-1);
fs=1/dt;

%f=(0:N-1)/N*fs;%frekvensakse
f_shift=( (-N/2) : 1 : (N/2-1) )*(fs/N);
P=fft(p);% fouriertransform av trykk p
P_shift=fftshift(P);

subplot(2,1,2);plot(f_shift,abs(P_shift));grid;%frekvensplot
%subplot(2,1,2);plot(f,abs(P));grid;%frekvensplot
title("Absolute pressure freq spectra")

%%
w=2*pi*f;%vinkelfrekvens
R=1; C=0.001;
Adm= 1/R - i*w*C;% admitans for paralellkobling av R og C
Z=1 ./Adm;% impedans

Q=Adm.*P;% flow i frekvensplan.  Fra Ohms lov P=Z*Q

q=real(ifft(Q));% flow i tidsplan

figure(2);clf;
subplot(2,1,1);plot(t,p);grid;
title("Pressure")
subplot(2,1,2);plot(t,q);grid;
title("Velocity")

