%% Patient18

addpath('Master/Prosjektoppgave/dataset/patient17/') 
load 20190117T145728_IQ_Sepsis-4min_traces;
name1 = '19.01.2019'; Ts1 = Ts; Tmean1 = Tmean; Tmin1 = Tmin; Tmax1 = Tmax; delay1 = -0.18; 
load 20190118T110902_IQ_Sepsis-4min_traces;
name2 = '20.01.2019'; Ts2 = Ts; Tmean2 = Tmean; Tmin2 = Tmin; Tmax2 = Tmax; delay2 = -0.18;
load 20190120T111256_IQ_Sepsis-4min_traces;
name3 = '21.01.2019'; Ts3 = Ts; Tmean3 = Tmean; Tmin3 = Tmin; Tmax3 = Tmax; delay3 = 0.07; 
load 20190123T105641_IQ_Sepsis-4min_traces;
name4 = '23.01.2019'; Ts4 = Ts; Tmean4 = Tmean; Tmin4 = Tmin; Tmax4 = Tmax; delay4 = -0.15; 
%% Register heart pulses on our own?
import register_heart_pulses.*

%t_pulses1 = register_heart_pulses(Ts1.ecg, Ts1.tED, 0.4);
%t_pulses2 = register_heart_pulses(Ts2.ecg, Ts2.tED, 0.4);
%t_pulses3 = register_heart_pulses(Ts3.ecg, Ts3.tED, 0.4);
%t_pulses4 = register_heart_pulses(Ts4.ecg, Ts4.tED, 0.4);

%% Else

t_pulses1 = Tmean1.tED;
t_pulses2 = Tmean2.tED;
t_pulses3 = Tmean3.tED;
t_pulses4 = Tmean4.tED;

%% Get compliance and resistance for all dates using fit

addpath('Master/Prosjektoppgave/helper_functions/') 
import calc_parameters.*

dataset.Ts=Ts1; dataset.Tmean=Tmean1; dataset.delay=delay1; dataset.t_pulses=t_pulses1;
[R1, C1, T1] = calc_parameters(dataset);
dataset.Ts=Ts2; dataset.Tmean=Tmean2; dataset.delay=delay2; dataset.t_pulses=t_pulses2;
[R2, C2, T2] = calc_parameters(dataset);
dataset.Ts=Ts3; dataset.Tmean=Tmean3; dataset.delay=delay3; dataset.t_pulses=t_pulses3;
[R3, C3, T3] = calc_parameters(dataset);
dataset.Ts=Ts4; dataset.Tmean=Tmean4; dataset.delay=delay4; dataset.t_pulses=t_pulses4;
[R4, C4, T4] = calc_parameters(dataset);

%% Interpolate to upsample and get a linear time axis using PCHIP (Piecewise Cubic Hermite Interpolating Polynomial)
interp_fac = 10;

T1_intp = linspace(T1(1), T1(end), length(T1)*interp_fac); 
T2_intp = linspace(T2(1), T2(end), length(T2)*interp_fac); 
T3_intp = linspace(T3(1), T3(end), length(T3)*interp_fac); 
T4_intp = linspace(T4(1), T4(end), length(T4)*interp_fac); 

R1_intp = interp1(T1,R1,T1_intp,'PCHIP'); C1_intp = interp1(T1,C1,T1_intp,'PCHIP');
R2_intp = interp1(T2,R2,T2_intp,'PCHIP'); C2_intp = interp1(T2,C2,T2_intp,'PCHIP');
R3_intp = interp1(T3,R3,T3_intp,'PCHIP'); C3_intp = interp1(T3,C3,T3_intp,'PCHIP');
R4_intp = interp1(T4,R4,T4_intp,'PCHIP'); C4_intp = interp1(T4,C4,T4_intp,'PCHIP');

N1=length(T1_intp);Tsamp=T1_intp(2)-T1_intp(1);fs=1/Tsamp;f1=(0:1/(N1-1):1)*fs;
N2=length(T2_intp);Tsamp=T2_intp(2)-T2_intp(1);fs=1/Tsamp;f2=(0:1/(N2-1):1)*fs;
N3=length(T3_intp);Tsamp=T3_intp(2)-T3_intp(1);fs=1/Tsamp;f3=(0:1/(N3-1):1)*fs;
N4=length(T4_intp);Tsamp=T4_intp(2)-T4_intp(1);fs=1/Tsamp;f4=(0:1/(N4-1):1)*fs;

%% Plot resistance

% Remove discontinuity between start/end
R1_intp = R1_intp - linspace(R1_intp(1), R1_intp(end), length(R1_intp)) + mean([R1_intp(1), R1_intp(end)]);
R2_intp = R2_intp - linspace(R2_intp(1), R2_intp(end), length(R2_intp)) + mean([R2_intp(1), R2_intp(end)]);
R3_intp = R3_intp - linspace(R3_intp(1), R3_intp(end), length(R3_intp)) + mean([R3_intp(1), R3_intp(end)]);
R4_intp = R4_intp - linspace(R4_intp(1), R4_intp(end), length(R4_intp)) + mean([R4_intp(1), R4_intp(end)]);

figure(10);clf;sgtitle('Patient 17 - Resistance');
subplot(4,1,1);plot(T1_intp,R1_intp);ylim([2000,7000]);grid();title(name1); set(gca,'XTickLabel',[]);
subplot(4,1,2);plot(T2_intp,R2_intp);ylim([2000,3200]);grid();title(name2); set(gca,'XTickLabel',[]);
subplot(4,1,3);plot(T3_intp,R3_intp);ylim([-20000,100000]);grid();title(name3); set(gca,'XTickLabel',[]);
subplot(4,1,4);plot(T4_intp,R4_intp);ylim([700,1200]);grid();title(name4); xlabel('Time [sec]')

%% Test DFT of resistance

% Remove (decrease) the DC component to reduce spectral bleeding from its
% ripples by subtracting the avg resistance and Normalize it and divide it 
% by DC component to get a similar relative scale for all of the spectra
R_dft1 = fft(R1_intp-mean(R1_intp))/N1;
R_dft2 = fft(R2_intp-mean(R2_intp))/N2;
R_dft3 = fft(R3_intp-mean(R3_intp))/N3; 
R_dft4 = fft(R4_intp-mean(R4_intp))/N4;

% Save data for workspace
R_dft1_patient17 = R_dft1;
R_dft2_patient17 = R_dft2;
R_dft3_patient17 = R_dft3;
R_dft4_patient17 = R_dft4;

figure(11);clf;sgtitle('Patient 17 - Resistance DFT');
subplot(4,2,1);plot(f1,abs(R_dft1));title(name1);xlim([0 0.1]);ylim([0 200]);grid();
subplot(4,2,3);plot(f2,abs(R_dft2));title(name2);xlim([0 0.1]);ylim([0 200]);grid();
subplot(4,2,5);plot(f3,abs(R_dft3));title(name3);xlim([0 0.1]);ylim([0 20000]);grid();
subplot(4,2,7);plot(f4,abs(R_dft4));title(name4);xlim([0 0.1]);ylim([0 200]);grid();

% Convert it to a comparable scale by converting values relative to their
% original DC conponent
R_dft1 = R_dft1/mean(R1_intp);
R_dft2 = R_dft2/mean(R2_intp);
R_dft3 = R_dft3/mean(R3_intp);
R_dft4 = R_dft4/mean(R4_intp);

% Save data for workspace
R_dft1_patient17 = R_dft1;
R_dft2_patient17 = R_dft2;
R_dft3_patient17 = R_dft3;
R_dft4_patient17 = R_dft4;

subplot(4,2,2);plot(f1,abs(R_dft1));title(name1);xlim([0 0.1]);ylim([0 3e-2]);grid();
subplot(4,2,4);plot(f2,abs(R_dft2));title(name2);xlim([0 0.1]);ylim([0 3e-2]);grid();
subplot(4,2,6);plot(f3,abs(R_dft3));title(name3);xlim([0 0.1]);ylim([0 3e-2]);grid();
subplot(4,2,8);plot(f4,abs(R_dft4));title(name4);xlim([0 0.1]);ylim([0 3e-2]);grid();

%% Crosscorrelate my method with actual measured data. We will recreate the 
% signal by using the lower frequency components to see if they match up
%figure(15)
%hold on
%plot(T4_intp,R4_intp);grid();title(name4)
%dft = fft(R4_intp)/N4;
%phi = angle(dft(18)); A = abs(dft(18))*2; % *2 because; 1 side. 
%phi2 = angle(dft(4)); A2 = abs(dft(4))*2; % *2 because; 1 side. 
%x = A*cos(2*pi*f4(18)*T4_intp + phi) + A2*cos(2*pi*f4(4)*T4_intp + phi2) + dft(1);
%%x = A2*cos(2*pi*f(6)*T_intp + phi2) + R_dft4(1)
%plot(T4_intp,x);ylim([0,2000]);grid();title(name4);grid();ylim([400,4000])
%hold off

%% Plot Capacitance

% Remove discontinuity between start/end
C1_intp = C1_intp - linspace(C1_intp(1), C1_intp(end), length(C1_intp)) + mean([C1_intp(1), C1_intp(end)]);
C2_intp = C2_intp - linspace(C2_intp(1), C2_intp(end), length(C2_intp)) + mean([C2_intp(1), C2_intp(end)]);
C3_intp = C3_intp - linspace(C3_intp(1), C3_intp(end), length(C3_intp)) + mean([C3_intp(1), C3_intp(end)]);
C4_intp = C4_intp - linspace(C4_intp(1), C4_intp(end), length(C4_intp)) + mean([C4_intp(1), C4_intp(end)]);

figure(12);clf;sgtitle('Patient 17 - Compliance');
subplot(4,1,1);plot(T1_intp,C1_intp);ylim([1.5e-4,2.5e-4]);grid();title(name1)
subplot(4,1,2);plot(T2_intp,C2_intp);ylim([2.8e-4,3.3e-4]);grid();title(name2)
subplot(4,1,3);plot(T3_intp,C3_intp);ylim([1e-4,3.3e-4]);grid();title(name3)
subplot(4,1,4);plot(T4_intp,C4_intp);ylim([1.8e-4,2.3e-4]);grid();title(name4)

%% Get DFT of Capacitance

% Remove (decrease) the DC component to reduce spectral bleeding from its
% ripples by subtracting the avg resistance and Normalize it and divide it 
% by DC component to get a similar relative scale for all of the spectra
C_dft1 = fft(C1_intp-mean(C1_intp))/N1;
C_dft2 = fft(C2_intp-mean(C2_intp))/N2;
C_dft3 = fft(C3_intp-mean(C3_intp))/N3; 
C_dft4 = fft(C4_intp-mean(C4_intp))/N4; 


figure(14);clf;sgtitle('Patient 18 - Compliance DFT');
subplot(4,2,1);plot(f1,abs(C_dft1));title(name1);xlim([0 0.1]);ylim([0 1.5e-5]);grid();
subplot(4,2,3);plot(f2,abs(C_dft2));title(name2);xlim([0 0.1]);ylim([0 1.5e-5]);grid();
subplot(4,2,5);plot(f3,abs(C_dft3));title(name3);xlim([0 0.1]);ylim([0 1.5e-5]);grid();
subplot(4,2,7);plot(f4,abs(C_dft4));title(name4);xlim([0 0.1]);ylim([0 1.5e-5]);grid();

% Convert it to a comparable scale by converting values relative to their
% original DC conponent
C_dft1 = C_dft1/mean(C1_intp);
C_dft2 = C_dft2/mean(C2_intp);
C_dft3 = C_dft3/mean(C3_intp);
C_dft4 = C_dft4/mean(C4_intp);

% Save data for workspace
C_dft1_patient17 = C_dft1;
C_dft2_patient17 = C_dft2;
C_dft3_patient17 = C_dft3;
C_dft4_patient17 = C_dft4;

subplot(4,2,2);plot(f1,abs(C_dft1));title(name1);xlim([0 0.1]);ylim([0 3e-2]);grid();
subplot(4,2,4);plot(f2,abs(C_dft2));title(name2);xlim([0 0.1]);ylim([0 3e-2]);grid();
subplot(4,2,6);plot(f3,abs(C_dft3));title(name3);xlim([0 0.1]);ylim([0 3e-2]);grid();
subplot(4,2,8);plot(f4,abs(C_dft4));title(name4);xlim([0 0.1]);ylim([0 3e-2]);grid();

