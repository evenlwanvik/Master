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
tIdx = find( Ts.t>tED(n) & Ts.t<tED(n+5)); % all samples within tED1 - tED2 window
t=Ts.t(tIdx);
v=Ts.velocity(tIdx);
p=Ts.ART(tIdx);

subplot(1,2,1);plot(t,v);title('Heart cycle - velocity')
subplot(1,2,2);plot(t,p);title('Heart cycle - pressure')

%% Plot velocity and pressure DFT for heart cycle

figure();clf;

N = length(t)
Ts = t(2)-t(1); fs = 1/Ts
f = (0:1/(N-1):1)*fs%
length(f)
length(p)
V = fft(v);
subplot(3,1,1);plot(f, abs(V));title('Heart cycle - velocity DFT')
xlim([0 fs/2]) % only show positive half of spectrum
P = fft(p);
subplot(3,1,2);plot(f, abs(P));title('Heart cycle - pressure DFT')
xlim([0 fs/2])

% Impedans for paralellkobling av R og C i Windkessel modell
Z = P/V;
R = Z(1) % zero freq component is the real value of Z
subplot(3,1,3);plot(f, abs(Z));title('Heart cycle - impedance DFT')
xlim([0 fs/2])

% find the first harmonic frequency of our signal
[x, i_Wv] = max(V(2:N-1)); Wv = f(i_Wv-1)
[x, i_Wp] = max(P(2:N-1)); Wp = f(i_Wp-1)
[x, i_Wz] = max(Z(2:N-1)); Wz = f(i_Wz-1)
1/Z(i_Wz)
Cw = (1/Wv)*(1/R-1/Z(i_Wz))
Cp = (1/Wp)*(1/R-1/Z(i_Wz))

%% Find a most likely capacitance using windkessel method

% parameters for 2 element
R = .95000 ; % systemic peripheral resistance i n (mmHg/cmˆ3/ sec )
C = 1.0666 ; % systemic arterial compliance i n (cmˆ3/mmHg)

%%Asumpltions
Tc= 60/72 ;% 72 beats per second
Ts =(2/5)*Tc ; % systole period
cycle =5; % number o f cardiac cycles for which WM i s analysed

% Modelling blood flow t o the aorta
syms ti q
I0 = solve( 90-int(q*(sin(pi*ti/Ts)),ti,0,Ts) , q );
I0 = subs(I0, '3.14', pi);
sine = @(t)sin(pi*t/Ts);
I = @(t)I0*sine(t).*(t<=Ts); % for one cycle

%% Plot all heart cycles to check if they all fit inside the tED windows

delay=0.07;
tED=Tmean.tED-delay;
plot_cycles(Ts, tED)

function plot_cycles(Ts, tED)
    figure(); clf();
    N = length(tED);
    for i=1:N
        idx = find( Ts.t>tED(i) & Ts.t<tED(i+1));
        subplot(1,2,1); plot(Ts.t(idx), Ts.velocity(idx));
        subplot(1,2,2); plot(Ts.t(idx), Ts.ART(idx));
        pause(1);
    end
end

