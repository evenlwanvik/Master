%% Patient18

addpath('Master/Prosjektoppgave/dataset/patient29/') 

load 20190820T095022_IQ_Sepsis-4min_traces;
name1 = '20.08.2019'; Ts1 = Ts; Tmean1 = Tmean; Tmin1 = Tmin; Tmax1 = Tmax; delay1=0.20;

load 20190821T094833_IQ_Sepsis-4min_traces;
name2 = '21.08.2019'; Ts2 = Ts; Tmean2 = Tmean; Tmin2 = Tmin; Tmax2 = Tmax; delay2=0.20;

load 20190821T164933_IQ_Sepsis-4min_traces;
name3 = '21.08.2019'; Ts3 = Ts; Tmean3 = Tmean; Tmin3 = Tmin; Tmax3 = Tmax; delay3=0.20;

load 20190823T083952_IQ_Sepsis-4min_traces;
name4 = '23.08.2019'; Ts4 = Ts; Tmean4 = Tmean; Tmin4 = Tmin; Tmax4 = Tmax; delay4=0.20;

load 20190823T134855_IQ_Sepsis-4min_traces;
name5 = '23.08.2019'; Ts5 = Ts; Tmean5 = Tmean; Tmin5 = Tmin; Tmax5 = Tmax; delay5=0.20;

load 20190826T103731_IQ_Sepsis-4min_traces;
name6 = '26.08.2019'; Ts6 = Ts; Tmean6 = Tmean; Tmin6 = Tmin; Tmax6 = Tmax; delay6=0.20;

load 20190827T120656_IQ_Sepsis-4min_traces;
name7 = '27.08.2019'; Ts7 = Ts; Tmean7 = Tmean; Tmin7 = Tmin; Tmax7 = Tmax; delay7=0.20; 


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
t_pulses5 = Tmean5.tED;
t_pulses6 = Tmean6.tED;
t_pulses7 = Tmean7.tED;

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
dataset.Ts=Ts5; dataset.Tmean=Tmean5; dataset.delay=delay5; dataset.t_pulses=t_pulses5;
[R5, C5, T5] = calc_parameters(dataset);
dataset.Ts=Ts6; dataset.Tmean=Tmean6; dataset.delay=delay6; dataset.t_pulses=t_pulses6;
[R6, C6, T6] = calc_parameters(dataset);
dataset.Ts=Ts7; dataset.Tmean=Tmean7; dataset.delay=delay7; dataset.t_pulses=t_pulses7;
[R7, C7, T7] = calc_parameters(dataset);




%% Interpolate to upsample and get a linear time axis using PCHIP (Piecewise Cubic Hermite Interpolating Polynomial)
interp_fac = 10;

T1_intp = linspace(T1(1), T1(end), length(T1)*interp_fac); 
T2_intp = linspace(T2(1), T2(end), length(T2)*interp_fac); 
T3_intp = linspace(T3(1), T3(end), length(T3)*interp_fac); 
T4_intp = linspace(T4(1), T4(end), length(T4)*interp_fac); 
T5_intp = linspace(T5(1), T5(end), length(T5)*interp_fac); 
T6_intp = linspace(T6(1), T6(end), length(T6)*interp_fac); 
T7_intp = linspace(T7(1), T7(end), length(T7)*interp_fac); 

R1_intp = interp1(T1,R1,T1_intp,'PCHIP'); C1_intp = interp1(T1,C1,T1_intp,'PCHIP');
R2_intp = interp1(T2,R2,T2_intp,'PCHIP'); C2_intp = interp1(T2,C2,T2_intp,'PCHIP');
R3_intp = interp1(T3,R3,T3_intp,'PCHIP'); C3_intp = interp1(T3,C3,T3_intp,'PCHIP');
R4_intp = interp1(T4,R4,T4_intp,'PCHIP'); C4_intp = interp1(T4,C4,T4_intp,'PCHIP');
R5_intp = interp1(T5,R5,T5_intp,'PCHIP'); C5_intp = interp1(T5,C5,T5_intp,'PCHIP');
R6_intp = interp1(T6,R6,T6_intp,'PCHIP'); C6_intp = interp1(T6,C6,T6_intp,'PCHIP');
R7_intp = interp1(T7,R7,T7_intp,'PCHIP'); C7_intp = interp1(T7,C7,T7_intp,'PCHIP');


N1=length(T1_intp);Tsamp=T1_intp(2)-T1_intp(1);fs=1/Tsamp;f1=(0:1/(N1-1):1)*fs;
N2=length(T2_intp);Tsamp=T2_intp(2)-T2_intp(1);fs=1/Tsamp;f2=(0:1/(N2-1):1)*fs;
N3=length(T3_intp);Tsamp=T3_intp(2)-T3_intp(1);fs=1/Tsamp;f3=(0:1/(N3-1):1)*fs;
N4=length(T4_intp);Tsamp=T4_intp(2)-T4_intp(1);fs=1/Tsamp;f4=(0:1/(N4-1):1)*fs;
N5=length(T5_intp);Tsamp=T5_intp(2)-T5_intp(1);fs=1/Tsamp;f5=(0:1/(N5-1):1)*fs;
N6=length(T6_intp);Tsamp=T6_intp(2)-T6_intp(1);fs=1/Tsamp;f6=(0:1/(N6-1):1)*fs;
N7=length(T7_intp);Tsamp=T7_intp(2)-T7_intp(1);fs=1/Tsamp;f7=(0:1/(N7-1):1)*fs;

%% Plot resistance

%figure(10);clf;sgtitle('Patient 18 - Resistance');
subplot(7,1,1);plot(T1_intp,R1_intp);%ylim([3000,6000]);grid();title(name1); set(gca,'XTickLabel',[]); %ylabel('SVR')
subplot(7,1,2);plot(T2_intp,R2_intp);%ylim([1500,3500]);grid();title(name2); set(gca,'XTickLabel',[]);
subplot(7,1,3);plot(T3_intp,R3_intp);%ylim([1000,3000]);grid();title(name3); set(gca,'XTickLabel',[]);
subplot(7,1,4);plot(T4_intp,R4_intp);%ylim([900,1400]);grid();title(name4); set(gca,'XTickLabel',[]);
subplot(7,1,5);plot(T5_intp,R5_intp);%ylim([900,1400]);grid();title(name5); set(gca,'XTickLabel',[]);
subplot(7,1,6);plot(T6_intp,R6_intp);%ylim([900,1400]);grid();title(name6); set(gca,'XTickLabel',[]);
subplot(7,1,7);plot(T7_intp,R7_intp);%ylim([900,1400]);grid();title(name7); xlabel('Time [s]')

%% Test relativity measure for Resistance
figure(70);clf;%sgtitle('Patient 18 - Resistance');
subplot(7,1,1);plot(T1_intp,R1_intp/mean(R1_intp));ylim([0.8,1.3]);grid();title(name1)
subplot(7,1,2);plot(T2_intp,R2_intp/mean(R2_intp));ylim([0.8,1.3]);grid();title(name2)
subplot(7,1,3);plot(T3_intp,R3_intp/mean(R3_intp));ylim([0.8,1.3]);grid();title(name3)
subplot(7,1,4);plot(T4_intp,R4_intp/mean(R4_intp));ylim([0.8,1.3]);grid();title(name4)
subplot(7,1,5);plot(T5_intp,R5_intp/mean(R5_intp));ylim([0.8,1.3]);grid();title(name5)
subplot(7,1,6);plot(T6_intp,R6_intp/mean(R6_intp));ylim([0.8,1.3]);grid();title(name6)
subplot(7,1,7);plot(T7_intp,R7_intp/mean(R7_intp));ylim([0.8,1.3]);grid();title(name7)


%% Plot relative Capacitance

figure(71);clf;%sgtitle('Patient 18 - Compliance');
subplot(7,1,1);plot(T1_intp,C1_intp/mean(C1_intp));ylim([0.7,1.5]);grid();title(name1)
subplot(7,1,2);plot(T2_intp,C2_intp/mean(C2_intp));ylim([0.7,1.5]);grid();title(name2)
subplot(7,1,3);plot(T3_intp,C3_intp/mean(C3_intp));ylim([0.7,1.5]);grid();title(name3)
subplot(7,1,4);plot(T4_intp,C4_intp/mean(C4_intp));ylim([0.7,1.5]);grid();title(name4)
subplot(7,1,5);plot(T5_intp,C5_intp/mean(C5_intp));ylim([0.7,1.5]);grid();title(name5)
subplot(7,1,6);plot(T6_intp,C6_intp/mean(C6_intp));ylim([0.7,1.5]);grid();title(name6)
subplot(7,1,7);plot(T7_intp,C7_intp/mean(C7_intp));ylim([0.7,1.5]);grid();title(name7)

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
R_dft5 = fft(R5_intp-mean(R5_intp))/N5/mean(R5_intp);
R_dft6 = fft(R6_intp-mean(R6_intp))/N6/mean(R6_intp);
R_dft7 = fft(R7_intp-mean(R7_intp))/N7/mean(R7_intp);

C_dft1 = fft(C1_intp-mean(C1_intp))/N1/mean(C1_intp);
C_dft2 = fft(C2_intp-mean(C2_intp))/N2/mean(C2_intp);
C_dft3 = fft(C3_intp-mean(C3_intp))/N3/mean(C3_intp);
C_dft4 = fft(C4_intp-mean(C4_intp))/N4/mean(C4_intp);
C_dft5 = fft(C5_intp-mean(C5_intp))/N5/mean(C5_intp);
C_dft6 = fft(C6_intp-mean(C6_intp))/N6/mean(C6_intp);
C_dft7 = fft(C7_intp-mean(C7_intp))/N7/mean(C7_intp);


% Extract elements from 2 (without DC component) to endIdx
R1_dft_subset = abs(R_dft1(idx));
R2_dft_subset = abs(R_dft2(idx));
R3_dft_subset = abs(R_dft3(idx));
R4_dft_subset = abs(R_dft4(idx));
R5_dft_subset = abs(R_dft5(idx));
R6_dft_subset = abs(R_dft6(idx));
R7_dft_subset = abs(R_dft7(idx));

% find average
R1_dft_subset_mean = mean(R1_dft_subset); 
R2_dft_subset_mean = mean(R2_dft_subset); 
R3_dft_subset_mean = mean(R3_dft_subset); 
R4_dft_subset_mean = mean(R4_dft_subset);
R5_dft_subset_mean = mean(R5_dft_subset);
R6_dft_subset_mean = mean(R6_dft_subset);
R7_dft_subset_mean = mean(R7_dft_subset);

R1_dft_subset_var = mad(R1_dft_subset);
R2_dft_subset_var = mad(R2_dft_subset);
R3_dft_subset_var = mad(R3_dft_subset);
R4_dft_subset_var = mad(R4_dft_subset);
R5_dft_subset_var = mad(R5_dft_subset);
R6_dft_subset_var = mad(R6_dft_subset);
R7_dft_subset_var = mad(R7_dft_subset);

R_mean = [R1_dft_subset_mean, R2_dft_subset_mean, R3_dft_subset_mean, R4_dft_subset_mean, R5_dft_subset_mean, R6_dft_subset_mean, R7_dft_subset_mean];
R_var = [R1_dft_subset_var, R2_dft_subset_var, R3_dft_subset_var, R4_dft_subset_var, R5_dft_subset_var, R6_dft_subset_var, R7_dft_subset_var];

plot(f_subset,R1_dft_subset, '-+', f_subset,R2_dft_subset, '-o' , f_subset,R3_dft_subset, '-*', f_subset,R4_dft_subset,'-s', f_subset,R5_dft_subset,'-d', f_subset,R6_dft_subset,'-^', f_subset,R7_dft_subset,'-v');
legend(name1,name2,name3,name4,name5,name6,name7); xlabel('Frequency [Hz]'); ylabel('(1)');

figure();
errorbar((1:7),R_mean, R_var, '-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
%title('Mean amplitude (~20-240sec)');grid();
xlim([0.75,7.25]);xticklabels({name1,name2,name3,name4,name5,name6,name7}); xlabel('Date of measurement'); ylabel('(1)');

%% Compliance


% Extract elements from 2 (without DC component) to endIdx
C1_dft_subset = abs(C_dft1(idx));
C2_dft_subset = abs(C_dft2(idx));
C3_dft_subset = abs(C_dft3(idx));
C4_dft_subset = abs(C_dft4(idx));
C5_dft_subset = abs(C_dft5(idx));
C6_dft_subset = abs(C_dft6(idx));
C7_dft_subset = abs(C_dft7(idx));

% find average
C1_dft_subset_mean = mean(C1_dft_subset); 
C2_dft_subset_mean = mean(C2_dft_subset); 
C3_dft_subset_mean = mean(C3_dft_subset); 
C4_dft_subset_mean = mean(C4_dft_subset);
C5_dft_subset_mean = mean(C5_dft_subset);
C6_dft_subset_mean = mean(C6_dft_subset);
C7_dft_subset_mean = mean(C7_dft_subset);

C1_dft_subset_var = mad(C1_dft_subset);
C2_dft_subset_var = mad(C2_dft_subset);
C3_dft_subset_var = mad(C3_dft_subset);
C4_dft_subset_var = mad(C4_dft_subset);
C5_dft_subset_var = mad(C5_dft_subset);
C6_dft_subset_var = mad(C6_dft_subset);
C7_dft_subset_var = mad(C7_dft_subset);

C_mean = [C1_dft_subset_mean, C2_dft_subset_mean, C3_dft_subset_mean, C4_dft_subset_mean, C5_dft_subset_mean, C6_dft_subset_mean, C7_dft_subset_mean];
C_var = [C1_dft_subset_var, C2_dft_subset_var, C3_dft_subset_var, C4_dft_subset_var, C5_dft_subset_var, C6_dft_subset_var, C7_dft_subset_var];
figure();
plot(f_subset,C1_dft_subset, '-+', f_subset,C2_dft_subset, '-o' , f_subset,C3_dft_subset, '-*', f_subset,C4_dft_subset,'-s', f_subset,C5_dft_subset,'-d', f_subset,C6_dft_subset,'-^', f_subset,C7_dft_subset,'-v');
legend(name1,name2,name3,name4,name5,name6,name7); xlabel('Frequency [Hz]'); ylabel('(1)');

figure();
errorbar((1:7),C_mean, C_var, '-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
%title('Mean amplitude (~20-240sec)');grid();
xlim([0.75,7.25]);xticklabels({name1,name2,name3,name4,name5,name6,name7}); xlabel('Date of measurement'); ylabel('(1)');

