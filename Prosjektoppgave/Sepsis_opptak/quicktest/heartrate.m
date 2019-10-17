%% Extract parameters over time for two-element Windkessel model.

addpath('Master/Prosjektoppgave/Sepsis_opptak/patient17/data') 
load 20190117T145728_IQ_Sepsis-4min_traces;
hr1 = Tmean.HR; t1 = Tmean.tED;
load 20190123T105641_IQ_Sepsis-4min_traces;
hr2 = Tmean.HR; t2 = Tmean.tED;

plot(t1,hr1,t2,hr2);
legend('17.11.2019', '23.01.2019');
grid();