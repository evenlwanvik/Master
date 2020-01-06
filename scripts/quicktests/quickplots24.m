% import functions
addpath('Master/Prosjektoppgave/helper_functions/') 
% add data
addpath('Master/Prosjektoppgave/dataset/patient24/') 
%load 20190319T090217_IQ_Sepsis-4min_traces; name = '19.03.2019'; delay=0.20;
%load 20190320T102140_IQ_Sepsis-4min_traces; name = '20.03.2019'; delay=0.20;
%load 20190321T111423_IQ_Sepsis-4min_traces; name = '21.03.2019'; delay=0.20;
load 20190323T161332_IQ_Sepsis-4min_traces; name = '23.03.2019'; delay=0.20; 

%% Find appropriate limit for registering heart pulses from ecg

figure();plot(Ts.t, Ts.ecg)

%% Register heart pulses on our own?
import register_heart_pulses.*

t_pulses = register_heart_pulses(Ts.ecg, Ts.t, 0.4);

%% test accuracy 
% by dividing the line from start to stop by a rato  and checking. 
% If it is higher than the diff between max and min we say that most likely
% the pulse registration is wrong.

dataset.Ts       = Ts; 
dataset.Tmean    = Tmean;
dataset.Tmax     = Tmax;
dataset.Tmin     = Tmin;
dataset.delay    = delay;
% We use Tmean.tED if every pulse is registered correctly
dataset.t_pulses = Tmean.tED;
%dataset.t_pulses = t_pulses;

import test_accuracy.*
[velocityMissRate, pressureMissRate] = test_accuracy(dataset, 8)

%% Plot all pulses to quickly investigate the waveforms
dataset.Ts       = Ts; 
dataset.Tmean    = Tmean;
dataset.delay    = delay;
% We use Tmean.tED if every pulse is registered correctly
dataset.t_pulses = Tmean.tED;
%dataset.t_pulses = t_pulses;

import plot_pulses.*
plot_pulses(dataset, 0.1)

%% Plot all dft's of pulse measurements
dataset.Ts       = Ts; 
dataset.Tmean    = Tmean;
dataset.delay    = delay;
% We use Tmean.tED if every pulse is registered correctly
dataset.t_pulses = Tmean.tED;
%dataset.t_pulses = t_pulses;

import plot_dfts.*
plot_dfts(dataset, 0.3)

%% plot full velocity and pressure
figure();
yyaxis left
plot(Ts.t,Ts.velocity);ylabel('Velocity [cm/s]');xlabel('t [s]')
yyaxis right
plot(Ts.t,Ts.ART);ylabel('Pressure [mmHg]')
xlim([59.4 61])

%% Test frequency
tED = t_pulses-delay;
for i=105:120
tIdx = find( Ts.t>tED(109) & Ts.t<tED(109+1)); % all samples within tED1 - tED2 window
t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx); N=length(t);
Tsamp = t(2)-t(1); fs = 1/Tsamp; f = (0:1/(N-1):1)*fs;
v_line = linspace(v(1), v(end), N); v = v - v_line'+mean(v_line);  
p_line = linspace(p(1), p(end), N); p = p - p_line'+mean(p_line);
V = fft(v); P = fft(p); Z = P./V;
Z(1)
end
%%
figure()
subplot(2,1,1);
yyaxis left
plot(f, abs(P))
yyaxis right
plot(f, abs(V));grid();
subplot(2,1,2);
plot(f, abs(Z));grid();
xlim([0, 15])
%% Heartrate

plot(Tmean1.tED,Tmean1.HR,Tmean2.tED,Tmean2.HR,Tmean3.tED,Tmean3.HR,Tmean4.tED,Tmean4.HR);
title('Heartrate before and after sepsis'); legend(name1,name2,name3,name4);
grid();ylim([240 90]);


%% Mean pressure

plot(Tmax1.tED,Tmax1.ART,Tmax2.tED,Tmax2.ART,Tmax3.tED,Tmax3.ART,Tmax4.tED,Tmax4.ART);
title('Mean pressure measurements'); legend('17.11.2019 1','17.11.2019 2','18.01.2019','20,01,2019','23.01.2019');
grid();ylim([0 160]);


%% Pressure

figure(); clf(); 
subplot(2,1,1); 
plot(Tmax1.tED,Tmax1.ART, Tmean1.tED,Tmean1.ART, Tmin1.tED,Tmin1.ART); ylim([0 160]); %xlim([5 20])
title('Pressure 17.11.2019'); legend('Max', 'Mean', 'Min'); grid();
subplot(2,1,2);  
plot(Tmax2.tED,Tmax2.ART, Tmean2.tED,Tmean2.ART, Tmin2.tED,Tmin2.ART); ylim([0 160]); %xlim([5 20])
title('Pressure 23.01.2019'); legend('Max', 'Mean', 'Min'); grid();

%% Filter the signal 
v = Ts.velocity; p = Ts.ART; t = Ts.t;
N = length(v); fs=1/(Ts.t(2)-Ts.t(1));
f = 0.5; fnorm = f/fs
d = designfilt('lowpassfir', ...
    'PassbandFrequency',0.7,'StopbandFrequency', 1, ...
    'PassbandRipple',1,'StopbandAttenuation',100, ...
    'DesignMethod','equiripple', 'SampleRate', fs);
figure(16);freqz(Num,1);
v_out=filtfilt(d,v);
p_out=filtfilt(d,p);
%%
del = 100; % number of elements to be deleted in beginning and end of arr
v_out = v_out(del+1:end-del);
p_out = p_out(del+1:end-del);
t = t(del+1:end-del);
figure(17); 
title('23.01.2019')
yyaxis left; plot(t,p_out);ylim([52, 70])
ylabel('pressure [mmHg]');
yyaxis right; plot(t,v_out);ylim([0.04, 0.082])%axis([0, 246, 0.040, 0.085])
ylabel('velocity [cm/s]');
xlabel('Time [sec]');
%% Get the complex envelope

figure(18);     
tED = Tmean.tED-delay;
tIdx = find( Ts.t>tED(50) & Ts.t<tED(50+1));
t = Ts.t(tIdx); p = Ts.ART(tIdx); v = Ts.velocity(tIdx);
P = abs(fft(p))/length(p);
fs = 1/(t(2)-t(1)); f = (0:1/(length(t)-1):1)*fs;
x = P(1);
for i=2:10
    x = x + 2*P(i)*exp(1j*2*pi*f(i).*(t+0.25));
end
plot(t,p,t,abs(x));%ylim([52, 70])
ylabel('pressure [mmHg]');
xlabel('Time [sec]')
legend('Measured pressure', 'Reconstructed pressure')
%%
yyaxis right; plot(Ts.t(tIdx),Ts.velocity(tIdx));%ylim([0.04, 0.082])%axis([0, 246, 0.040, 0.085])
ylabel('velocity [cm/s]');
xlabel('Time [sec]');