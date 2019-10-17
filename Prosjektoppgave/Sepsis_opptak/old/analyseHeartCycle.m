%% Analyse av pasient nr 17
% We have been given two snippets of measurements on different dates: 
% * 17/01/2019 (1) and 23/01/2019

load 20190117T145728_IQ_Sepsis-4min_traces;
Tmax1 = Tmax; Tmean1 = Tmean; Tmin1 = Tmin; Ts1 = Ts;
load 20190123T105641_IQ_Sepsis-4min_traces;
Tmax2 = Tmax; Tmean2 = Tmean; Tmin2 = Tmin; Ts2 = Ts;


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
tED=Tmean.tED-delay;
Ncycles=length(tED);

% Plot velocity and pressure for heart cycle
n=10;
tIdx = find( Ts.t>tED(n) & Ts.t<tED(n+1));
N = length(tIdx)
t=Ts.t(tIdx);
v=Ts.velocity(tIdx);
p=Ts.ART(tIdx);
subplot(2,2,1);plot(t,v);title('Heart cycle - velocity')
subplot(2,2,3);plot(t,p);title('Heart cycle - pressure')

%% Plot velocity and pressure DFT for heart cycle

figure();clf;

V = fft(v);
Vdb = mag2db(abs(V))
subplot(3,1,1);semilogx(Vdb);title('Heart cycle - velocity DFT')
P = fft(p);
subplot(3,1,2);semilogx(abs(P));title('Heart cycle - pressure DFT')

% Impedans for paralellkobling av R og C i Windkessel modell
Z = P/V;
subplot(3,1,3);semilogx(abs(Z));title('Heart cycle - impedance DFT')

%% Calculate resistance (DC component)