%% Winkessel demo
% Endre motstand R og kondensator C, og observer endring i flow kurveform
% Regn ut pulsativ index  (PI) for p og q
% Hvordan endrer PI seg med varierende R og C?

% 2019.08.19 . Hans Torp

%load 20190117T145728_IQ_Sepsis-4min_traces;
%date = '17.01.2019'
load 20190123T105641_IQ_Sepsis-4min_traces;
date = '23.01.2019'
name = append('Sepsis 4min_traces ', date)
figure();clf;sgtitle(name);

% extract a single heart cycle
delay = 0.07;
tED = Tmean.tED-delay;
cycleStart=33; % choose heart cycle (~300 possibilities)
nCycles=1; % choose how many heart cycles to use
tIdx = find( Ts.t>tED(cycleStart) & Ts.t<tED(cycleStart+nCycles)); % all samples within tED1 - tED2 window
t=Ts.t(tIdx);
v=Ts.velocity(tIdx);
p=Ts.ART(tIdx);

subplot(2,1,1);plot(t,p);grid;
title("Pressure")

N=length(t);
dt=(t(end)-t(1))/(N-1);
fs=1/dt;

f=(0:N-1)/N*fs;%frekvensakse
P=fft(p);% fouriertransform av trykk p
V=fft(v);% fouriertransform av trykk p
Z=P./V;

subplot(2,1,2);plot(f,abs(P));grid;%frekvensplot
xlim([0,10])
title("Absolute pressure freq spectra")

%%
w=2*pi*f;%vinkelfrekvens
R=1000; C=1.18e-4;
Adm= 1/R - i*w*C;% admitans for paralellkobling av R og C
Z_WK=1./Adm;% impedans

Q=Adm.*P;% flow i frekvensplan.  Fra Ohms lov P=Z*Q

q=real(ifft(Q));% flow i tidsplan

figure();clf;
subplot(3,1,1);plot(t,p);grid;
title("Pressure")
subplot(3,1,2);plot(t,v);grid;
title("Velocity")
subplot(3,1,3);plot(t,q);grid;
title("Velocity")

