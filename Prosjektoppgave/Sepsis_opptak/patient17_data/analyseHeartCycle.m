%% Analyse av pasient nr 17
% We have been given two snippets of measurements on different dates: 
% * 17/01/2019 (1) and 23/01/2019

load 20190117T145728_IQ_Sepsis-4min_traces;
Tmax1 = Tmax; Tmean1 = Tmean; Tmin1 = Tmin; Ts1 = Ts;
load 20190123T105641_IQ_Sepsis-4min_traces;
Tmax2 = Tmax; Tmean2 = Tmean; Tmin2 = Tmin; Ts2 = Ts;
T = Ts2.t(2)-Ts2.t(1); fs = 1/T;

%% Normal plot of both datasets

figure();clf;
subplot(2,2,1);plot(Ts1.t, Ts1.velocity);grid;
title("Raw velocity measurement - 20190117")
subplot(2,2,2);plot(Ts1.t, Ts1.ART);grid;
title("Raw pressure measurement - 20190117")
subplot(2,2,3);plot(Ts2.t, Ts2.velocity);grid;
title("Raw velocity measurement - 20190123")
subplot(2,2,4);plot(Ts2.t, Ts2.ART);grid;
title("Raw pressure measurement - 20190123")

%% pick one heart cycle
figure();clf;

delay=0.07;
tED=Tmean2.tED-delay;
Ncycles=length(tED);

% Plot velocity and pressure for heart cycle
n=10;
tIdx = find( Ts2.t>tED(n) & Ts2.t<tED(n+3)); % all samples within tED1 - tED2 window
t=Ts2.t(tIdx); 
v=Ts2.velocity(tIdx);
p=Ts2.ART(tIdx);

subplot(1,2,1);plot(t,v);title('Heart cycle - velocity')
subplot(1,2,2);plot(t,p);title('Heart cycle - pressure')

%% Extract parameters over time for two-element Windkessel model.

figure();clf;

N = length(t)
f = (0:1/(N-1):1)*fs;
V = fft(v);
subplot(3,1,1);plot(f, mag2db(abs(V)));title('Heart cycle - velocity DFT')
xlim([0 10]) % only show positive half of spectrum
P = fft(p);
subplot(3,1,2);plot(f, mag2db(abs(P)));title('Heart cycle - pressure DFT')
xlim([0 10])

size(P)
size(V.')

% Impedans for paralellkobling av R og C i Windkessel modell
Z = P./V;
%Z = rdivide(P,V)
R = Z(1) % zero freq component is the real value of Z
subplot(3,1,3);plot(f, mag2db(abs(Z)));title('Heart cycle - impedance DFT')
xlim([0 10])

% find the first harmonic frequency of our signal
[x, i_Fv] = max(abs(V(2:N-1))); Fv = f(i_Fv-1)
[x, i_Fp] = max(abs(P(2:N-1))); Fp = f(i_Fp-1)
[x, i_Fz] = max(abs(Z(2:N-1))); Fz = f(i_Fz-1)

%Cz = abs( Z(i_Fp)/(1i*Fp) )
Cz = abs( (1/R-1/Z(i_Fp))/i_Fz )
Adm = 1/R + i*Fp*0.01;
Q=Adm.*P;
q=abs(real(ifft(Q)));% flow i tidsplan
%subplot(1,3,1);plot(t,p);title('Heart cycle - pressure')
%subplot(1,3,2);plot(t,v);title('Heart cycle - velocity')
%subplot(1,3,3);plot(t,q);title('Velocity remade')


