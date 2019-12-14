%% Patient16
% import functions
addpath('Master/Prosjektoppgave/helper_functions/') 

% What dataset you wanna play with?
addpath('Master/Prosjektoppgave/dataset/patient16/') 
load 20181213T125033_IQ_Sepsis_traces; name = '13.12.2018'; 
%load 20181213T205123_IQ_Sepsis-4min_traces; name = '13.12.2018'; 
%load 20181214T090138_IQ_Sepsis-4min_traces; name = '14.12.2018'; 
%load 20181217T085556_IQ_Sepsis_traces; name = '17.12.2018';
delay = -0.15;

%% Find appropriate limit for registering heart pulses from ecg

figure();plot(Ts.t, Ts.ecg)

%% Register heart pulses on our own?
import register_heart_pulses.*

t_pulses = register_heart_pulses(Ts.ecg, Ts.t, 0.2);

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
dataset.t_pulses = Tmean.t;
dataset.t_pulses = t_pulses;

import test_accuracy.*
[velocityMissRate, pressureMissRate] = test_accuracy(dataset, 10)

%% Plot all pulses to quickly investigate the waveforms
dataset.Ts       = Ts; 
dataset.Tmean    = Tmean;
dataset.delay    = delay;
% We use Tmean.tED if every pulse is registered correctly
%dataset.t_pulses = Tmean.t;
dataset.t_pulses = t_pulses;

import plot_pulses.*
plot_pulses(dataset, 0.1)

%% Plot all dft's of pulse measurements
dataset.Ts       = Ts; 
dataset.Tmean    = Tmean;
dataset.delay    = delay;
% We use Tmean.tED if every pulse is registered correctly
dataset.t_pulses = Tmean.t;
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


%% Plot full
figure(3)
subplot(3,1,1);plot(Ts1.t, Ts1.velocity);title('Velocity 17.11.2019 1')
subplot(3,1,2);plot(Ts4.t, Ts4.velocity);grid;title('Velocity 20,01,2019')

figure(4)
%subplot(2,1,1);plot(Ts5.t, Ts5.velocity)
%subplot(2,1,2);plot(Ts5.t, Ts5.ART)