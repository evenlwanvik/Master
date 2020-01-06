addpath('Master/Prosjektoppgave/Sepsis_opptak/patient17/data') 
load 20190117T145728_IQ_Sepsis-4min_traces;
Ts1 = Ts; Tmean1 = Tmean; Tmin1 = Tmin; Tmax1 = Tmax; 
load 20190117T215434_IQ_Sepsis-4min_traces;
Ts2 = Ts; Tmean2 = Tmean; Tmin2 = Tmin; Tmax2 = Tmax; 
load 20190118T110902_IQ_Sepsis-4min_traces;
Ts3 = Ts; Tmean3 = Tmean; Tmin3 = Tmin; Tmax3 = Tmax; 
load 20190120T111256_IQ_Sepsis-4min_traces;
Ts4 = Ts; Tmean4 = Tmean; Tmin4 = Tmin; Tmax4 = Tmax; 
load 20190123T105641_IQ_Sepsis-4min_traces;
Ts5 = Ts; Tmean5 = Tmean; Tmin5 = Tmin; Tmax5 = Tmax; 

%% Heartrate
figure();
plot(Tmean1.tED,Tmean1.HR,Tmean2.tED,Tmean2.HR,Tmean3.tED,Tmean3.HR,Tmean4.tED,Tmean4.HR,Tmean5.tED,Tmean5.HR);
title('Heartrate before and after sepsis'); legend('17.11.2019 1','17.11.2019 2','18.01.2019','20,01,2019','23.01.2019');
grid();ylim([0 160]);

%% Mean pressure

plot(Tmax1.tED,Tmax1.ART,Tmax2.tED,Tmax2.ART,Tmax3.tED,Tmax3.ART,Tmax4.tED,Tmax4.ART,Tmax5.tED,Tmax5.ART);
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
figure()
plot(Tmean4.tED, Tmean4.velocity, '*', Ts4.t, Ts4.velocity )
figure()
plot(Tmean2.tED, Tmean2.velocity, '*', Ts2.t, Ts2.velocity )