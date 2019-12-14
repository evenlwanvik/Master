%% Patient 24

addpath('Master/Prosjektoppgave/dataset/patient24/') 
load 20190319T090217_IQ_Sepsis-4min_traces;
name1 = '19.03.2019'; Ts1 = Ts; Tmean1 = Tmean; Tmin1 = Tmin; Tmax1 = Tmax; delay1=0.20;
load 20190320T102140_IQ_Sepsis-4min_traces;
name2 = '20.03.2019'; Ts2 = Ts; Tmean2 = Tmean; Tmin2 = Tmin; Tmax2 = Tmax; delay2=0.20;
load 20190321T111423_IQ_Sepsis-4min_traces;
name3 = '21.03.2019'; Ts3 = Ts; Tmean3 = Tmean; Tmin3 = Tmin; Tmax3 = Tmax; delay3=0.20;
load 20190323T161332_IQ_Sepsis-4min_traces;
name4 = '23.03.2019'; Ts4 = Ts; Tmean4 = Tmean; Tmin4 = Tmin; Tmax4 = Tmax; delay4=0.20; 

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

%figure(10);clf;sgtitle('Patient 18 - Resistance');
subplot(4,1,1);plot(T1_intp,R1_intp);%ylim([3000,6000]);grid();title(name1); set(gca,'XTickLabel',[]); %ylabel('SVR')
subplot(4,1,2);plot(T2_intp,R2_intp);%ylim([1500,3500]);grid();title(name2); set(gca,'XTickLabel',[]);
subplot(4,1,3);plot(T3_intp,R3_intp);%ylim([1000,3000]);grid();title(name3); set(gca,'XTickLabel',[]);
subplot(4,1,4);plot(T4_intp,R4_intp);%ylim([900,1400]);grid();title(name4); set(gca,'XTickLabel',[]);

%% Test relativity measure for Resistance
figure(70);clf;%sgtitle('Patient 18 - Resistance');
subplot(4,1,1);plot(T1_intp,R1_intp/mean(R1_intp));ylim([0.8,1.2]);grid();title(name1)
subplot(4,1,2);plot(T2_intp,R2_intp/mean(R2_intp));ylim([0.8,1.2]);grid();title(name2)
subplot(4,1,3);plot(T3_intp,R3_intp/mean(R3_intp));ylim([0.8,1.2]);grid();title(name3)
subplot(4,1,4);plot(T4_intp,R4_intp/mean(R4_intp));ylim([0.8,1.2]);grid();title(name4)

%% Plot relative Capacitance

figure(71);clf;%sgtitle('Patient 18 - Compliance');
subplot(4,1,1);plot(T1_intp,C1_intp/mean(C1_intp));ylim([0.8,1.3]);grid();title(name1)
subplot(4,1,2);plot(T2_intp,C2_intp/mean(C2_intp));ylim([0.8,1.3]);grid();title(name2)
subplot(4,1,3);plot(T3_intp,C3_intp/mean(C3_intp));ylim([0.8,1.3]);grid();title(name3)
subplot(4,1,4);plot(T4_intp,C4_intp/mean(C4_intp));ylim([0.8,1.3]);grid();title(name4)

%% Plot lower frequencies


% we wants oscillations from 20 sec to 2min, i.e. 0-0.05  Hz
N1=length(T1_intp);Tsamp=T1_intp(2)-T1_intp(1);fs=1/Tsamp;f1=(0:1/(N1-1):1)*fs;
startIdx = 3; endIdx = round(0.04/(fs/N1))+1; idx = startIdx:endIdx;

f_subset = f1(idx);
figure(80);clf;

R_dft1 = fft(R1_intp-mean(R1_intp))/N1/mean(R1_intp);
R_dft2 = fft(R2_intp-mean(R2_intp))/N2/mean(R2_intp);
R_dft3 = fft(R3_intp-mean(R3_intp))/N3/mean(R3_intp);
R_dft4 = fft(R4_intp-mean(R4_intp))/N4/mean(R4_intp);

C_dft1 = fft(C1_intp-mean(C1_intp))/N1/mean(C1_intp);
C_dft2 = fft(C2_intp-mean(C2_intp))/N2/mean(C2_intp);
C_dft3 = fft(C3_intp-mean(C3_intp))/N3/mean(C3_intp);
C_dft4 = fft(C4_intp-mean(C4_intp))/N4/mean(C4_intp);

% Extract elements from 2 (without DC component) to endIdx
R1_dft_subset = abs(R_dft1(idx));
R2_dft_subset = abs(R_dft2(idx));
R3_dft_subset = abs(R_dft3(idx));
R4_dft_subset = abs(R_dft4(idx));

% find average
R1_dft_subset_mean = mean(R1_dft_subset); 
R2_dft_subset_mean = mean(R2_dft_subset); 
R3_dft_subset_mean = mean(R3_dft_subset); 
R4_dft_subset_mean = mean(R4_dft_subset);

R1_dft_subset_var = mad(R1_dft_subset);
R2_dft_subset_var = mad(R2_dft_subset);
R3_dft_subset_var = mad(R3_dft_subset);
R4_dft_subset_var = mad(R4_dft_subset);

R_mean = [R1_dft_subset_mean, R2_dft_subset_mean, R3_dft_subset_mean, R4_dft_subset_mean];
R_var = [R1_dft_subset_var, R2_dft_subset_var, R3_dft_subset_var, R4_dft_subset_var];

plot(f_subset,R1_dft_subset, '-+', f_subset,R2_dft_subset, '-o' , f_subset,R3_dft_subset, '-*', f_subset,R4_dft_subset);
legend(name1,name2,name3,name4); xlabel('Frequency [Hz]'); ylabel('(1)');

figure();
errorbar((1:4),R_mean, R_var, '-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
%title('Mean amplitude (~20-240sec)');grid();
xlim([0.75,4.25]);xticklabels({name1,'',name2,'',name3,'',name4}); xlabel('Date of measurement'); ylabel('(1)');

%% Compliance


% Extract elements from 2 (without DC component) to endIdx
C1_dft_subset = abs(C_dft1(idx));
C2_dft_subset = abs(C_dft2(idx));
C3_dft_subset = abs(C_dft3(idx));
C4_dft_subset = abs(C_dft4(idx));

% find average
C1_dft_subset_mean = mean(C1_dft_subset); 
C2_dft_subset_mean = mean(C2_dft_subset); 
C3_dft_subset_mean = mean(C3_dft_subset); 
C4_dft_subset_mean = mean(C4_dft_subset);

C1_dft_subset_var = mad(C1_dft_subset);
C2_dft_subset_var = mad(C2_dft_subset);
C3_dft_subset_var = mad(C3_dft_subset);
C4_dft_subset_var = mad(C4_dft_subset);

C_mean = [C1_dft_subset_mean, C2_dft_subset_mean, C3_dft_subset_mean, C4_dft_subset_mean];
C_var = [C1_dft_subset_var, C2_dft_subset_var, C3_dft_subset_var, C4_dft_subset_var];
figure();
plot(f_subset,C1_dft_subset, '-+', f_subset,C2_dft_subset, '-o' , f_subset,C3_dft_subset, '-*', f_subset,C4_dft_subset,'-s');
legend(name1,name2,name3,name4,name5,name6,name7); xlabel('Frequency [Hz]'); ylabel('(1)');

figure();
errorbar((1:4),C_mean, C_var, '-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
%title('Mean amplitude (~20-240sec)');grid();
xlim([0.75,4.25]);xticklabels({name1,'',name2,'',name3,'',name4}); xlabel('Date of measurement'); ylabel('(1)');

