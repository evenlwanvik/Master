%% Patient16
addpath('Master/Prosjektoppgave/helper_functions/') 
addpath('Master/Prosjektoppgave/dataset/patient16/') 
load 20181213T125033_IQ_Sepsis_traces;
name1 = '13.12.2018'; Ts1 = Ts; Tmean1 = Tmean; Tmin1 = Tmin; Tmax1 = Tmax; 
load 20181213T205123_IQ_Sepsis-4min_traces;
name2 = '13.12.2018'; Ts2 = Ts; Tmean2 = Tmean; Tmin2 = Tmin; Tmax2 = Tmax; 
load 20181214T090138_IQ_Sepsis-4min_traces;
name3 = '14.12.2018'; Ts3 = Ts; Tmean3 = Tmean; Tmin3 = Tmin; Tmax3 = Tmax; 
load 20181217T085556_IQ_Sepsis_traces;
name4 = '17.12.2018'; Ts4 = Ts; Tmean4 = Tmean; Tmin4 = Tmin; Tmax4 = Tmax; 
delay = -0.15;

%% Register heart pulses on our own?
import register_heart_pulses.*

t_pulses1 = register_heart_pulses(Ts1.ecg, Ts1.t, 0.2);
t_pulses2 = register_heart_pulses(Ts2.ecg, Ts2.t, 0.2);
t_pulses3 = register_heart_pulses(Ts3.ecg, Ts3.t, 0.2);
t_pulses4 = register_heart_pulses(Ts4.ecg, Ts4.t, 0.2);

%% Else

%t_pulses1 = Tmean1.t;
%t_pulses2 = Tmean2.t;
%t_pulses3 = Tmean3.t;
%t_pulses4 = Tmean4.t;

%% Get compliance and resistance for all dates using fit

import calc_parameters.*

dataset.Ts=Ts1; dataset.Tmean=Tmean1; dataset.delay=delay; dataset.t_pulses=t_pulses1;
[R1, C1, T1] = calc_parameters(dataset);
dataset.Ts=Ts2; dataset.Tmean=Tmean2; dataset.delay=delay; dataset.t_pulses=t_pulses2;
[R2, C2, T2] = calc_parameters(dataset);
dataset.Ts=Ts3; dataset.Tmean=Tmean3; dataset.delay=delay; dataset.t_pulses=t_pulses3;
[R3, C3, T3] = calc_parameters(dataset);
dataset.Ts=Ts4; dataset.Tmean=Tmean4; dataset.delay=delay; dataset.t_pulses=t_pulses4;
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


%% Plot resistance

figure(10);clf;sgtitle('Patient 16 - Resistance using first harmonic for 1-component WK');
subplot(4,1,1);plot(T1_intp,R1_intp);ylim([0,10000]);grid();title(name1)
subplot(4,1,2);plot(T2_intp,R2_intp);ylim([0,10000]);grid();title(name2)
subplot(4,1,3);plot(T3_intp,R3_intp);ylim([1300,2000]);grid();title(name3)
subplot(4,1,4);plot(T4_intp,R4_intp);ylim([1000,4000]);grid();title(name4)

%% Test DFT of resistance

N1=length(T1_intp);Tsamp=T1_intp(2)-T1_intp(1);fs=1/Tsamp;f1=(0:1/(N1-1):1)*fs;
N2=length(T2_intp);Tsamp=T2_intp(2)-T2_intp(1);fs=1/Tsamp;f2=(0:1/(N2-1):1)*fs;
N3=length(T3_intp);Tsamp=T3_intp(2)-T3_intp(1);fs=1/Tsamp;f3=(0:1/(N3-1):1)*fs;
N4=length(T4_intp);Tsamp=T4_intp(2)-T4_intp(1);fs=1/Tsamp;f4=(0:1/(N4-1):1)*fs;

% Normalize it and divide it by DC component to get a similar relative
% scale for all of the spectra
R_dft1 = fft(R1_intp)/N1; %R_dft1 = R_dft1/R_dft1(1);
R_dft2 = fft(R2_intp)/N2; %R_dft2 = R_dft2/R_dft2(1);
R_dft3 = fft(R3_intp)/N3; %R_dft3 = R_dft3/R_dft3(1);
R_dft4 = fft(R4_intp)/N4; %R_dft4 = R_dft4/R_dft4(1);

% Save data for workspace
R_dft1_patient16 = R_dft1;
R_dft2_patient16 = R_dft2;
R_dft3_patient16 = R_dft3;
R_dft4_patient16 = R_dft4;

figure(11);clf;sgtitle('Patient 16 - Resistance DFT first harmonic for 1-component WK');
subplot(4,2,1);plot(f1,abs(R_dft1));title(name1);xlim([0 0.1]);ylim([0 10000]);grid();
subplot(4,2,3);plot(f2,abs(R_dft2));title(name2);xlim([0 0.1]);ylim([0 3000]);grid();
subplot(4,2,5);plot(f3,abs(R_dft3));title(name3);xlim([0 0.1]);ylim([0 3000]);grid();
subplot(4,2,7);plot(f4,abs(R_dft4));title(name4);xlim([0 0.1]);ylim([0 3000]);grid();

R_dft1 = R_dft1/R_dft1(1);
R_dft2 = R_dft2/R_dft2(1);
R_dft3 = R_dft3/R_dft3(1);
R_dft4 = R_dft4/R_dft4(1);

% Save data for workspace
R_dft1_patient16 = R_dft1;
R_dft2_patient16 = R_dft2;
R_dft3_patient16 = R_dft3;
R_dft4_patient16 = R_dft4;

subplot(4,2,2);plot(f1,abs(R_dft1));title(name1);xlim([0 0.1]);ylim([0 1]);grid();
subplot(4,2,4);plot(f2,abs(R_dft2));title(name2);xlim([0 0.1]);ylim([0 0.3]);grid();
subplot(4,2,6);plot(f3,abs(R_dft3));title(name3);xlim([0 0.1]);ylim([0 0.3]);grid();
subplot(4,2,8);plot(f4,abs(R_dft4));title(name4);xlim([0 0.1]);ylim([0 0.3]);grid();

%subplot(4,2,2);plot(f,mag2db(abs(R_dft1)));xlim([0,fs/2]);title(name1);xlim([0 3]);ylim([-60 -20]);xlim([0 0.5])
%subplot(4,2,4);plot(f,mag2db(abs(R_dft2)));xlim([0,fs/2]);title(name2);xlim([0 3]);ylim([-60 -20]);xlim([0 0.5])
%subplot(4,2,6);plot(f,mag2db(abs(R_dft3)));xlim([0,fs/2]);title(name3);xlim([0 3]);ylim([-60 -20]);xlim([0 0.5])
%subplot(4,2,8);plot(f,mag2db(abs(R_dft4)));xlim([0,fs/2]);title(name4);xlim([0 3]);ylim([-60 -20]);xlim([0 0.5])

%%
figure(15)
hold on
plot(T4_intp,R4_intp);grid();title(name4)
dft = fft(R4_intp)/N4;
phi = angle(dft(18)); A = abs(dft(18))*2; % *2 because; 1 side. 
phi2 = angle(dft(4)); A2 = abs(dft(4))*2; % *2 because; 1 side. 
x = A*cos(2*pi*f4(18)*T4_intp + phi) + A2*cos(2*pi*f4(4)*T4_intp + phi2) + dft(1);
%x = A2*cos(2*pi*f(6)*T_intp + phi2) + R_dft4(1)
plot(T4_intp,x);ylim([0,2000]);grid();title(name4);grid();ylim([400,4000])
hold off

%% db scale

%figure(11);clf;sgtitle('Compliance DFT first harmonic for 1-component WK');
%subplot(4,1,1);plot(f1,mag2db(abs(C_dft1)));xlim([0,fs/2]);title(name1);%ylim([-100,100])
%subplot(4,1,2);plot(f2,mag2db(abs(C_dft2)));xlim([0,fs/2]);title(name2);%ylim([-100,100])
%subplot(4,1,3);plot(f3,mag2db(abs(C_dft3)));xlim([0,fs/2]);title(name3);%ylim([-100,100])
%subplot(4,1,4);plot(f4,mag2db(abs(C_dft4)));xlim([0,fs/2]);title(name4);%ylim([-100,100])

%% Plot compliance

figure(12);clf;sgtitle('Patient 16 - Compliance using first harmonic for 1-component WK');
C_line = linspace(C1_intp(1), C1_intp(end), N1); C1_fin=C1_intp-C_line+mean([C1_intp(1) C1_intp(end)]);
subplot(4,1,1);plot(T1_intp,C1_fin);grid();title(name1);ylim([3.3e-4,5e-4])

C_line = linspace(C2_intp(1), C2_intp(end), N2); C2_fin=C2_intp-C_line+mean([C2_intp(1) C2_intp(end)]);
subplot(4,1,2);plot(T2_intp,C2_fin);grid();title(name2);ylim([2.8e-4,5e-4])

C_line = linspace(C3_intp(1), C3_intp(end), N3); C3_fin=C3_intp-C_line+mean([C3_intp(1) C3_intp(end)]);
subplot(4,1,3);plot(T3_intp,C3_fin);grid();title(name3);ylim([4.5e-4,6.5e-4])

C_line = linspace(C4_intp(1), C4_intp(end), N4); C4_fin=C4_intp-C_line+mean([C4_intp(1) C4_intp(end)]);
subplot(4,1,4);plot(T4_intp,C4_fin);grid();title(name4);ylim([3e-4,6e-4])

%% Test DFT of compliance

% Normalize it and divide it by DC component to get a similar relative
% scale for all of the spectra
C_dft1 = fft(C1_fin)/N1; %C_dft1 = C_dft1/C_dft1(1);
C_dft2 = fft(C2_fin)/N2; %R_dft2 = R_dft2/R_dft2(1);
C_dft3 = fft(C3_fin)/N3; %R_dft3 = R_dft3/R_dft3(1);
C_dft4 = fft(C4_fin)/N4; %R_dft4 = R_dft4/R_dft4(1);

figure(14);clf;sgtitle('Patient 16 - Compliance DFT first harmonic for 1-component WK');
subplot(4,2,1);plot(f1,abs(C_dft1));xlim([0,fs/2]);title(name1);xlim([0 0.1]);ylim([0 2e-5])
subplot(4,2,3);plot(f2,abs(C_dft2));xlim([0,fs/2]);title(name2);xlim([0 0.1]);ylim([0 2e-5])
subplot(4,2,5);plot(f3,abs(C_dft3));xlim([0,fs/2]);title(name3);xlim([0 0.1]);ylim([0 2e-5])
subplot(4,2,7);plot(f4,abs(C_dft4));xlim([0,fs/2]);title(name4);xlim([0 0.1]);ylim([0 2e-5])

C_dft1 = C_dft1/C_dft1(1);
C_dft2 = C_dft2/C_dft2(1);
C_dft3 = C_dft3/C_dft3(1);
C_dft4 = C_dft4/C_dft4(1);

% Save data for workspace
C_dft1_patient16 = C_dft1;
C_dft2_patient16 = C_dft2;
C_dft3_patient16 = C_dft3;
C_dft4_patient16 = C_dft4;

subplot(4,2,2);plot(f1,abs(C_dft1));xlim([0,fs/2]);title(name1);xlim([0 0.1]);ylim([0 0.05])
subplot(4,2,4);plot(f2,abs(C_dft2));xlim([0,fs/2]);title(name2);xlim([0 0.1]);ylim([0 0.05])
subplot(4,2,6);plot(f3,abs(C_dft3));xlim([0,fs/2]);title(name3);xlim([0 0.1]);ylim([0 0.05])
subplot(4,2,8);plot(f4,abs(C_dft4));xlim([0,fs/2]);title(name4);xlim([0 0.1]);ylim([0 0.05])
