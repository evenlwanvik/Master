% Run "week_44/compliance_patient##.m first"

% we wants oscillations from 20 sec to 2min, i.e. 0-0.05  Hz
N1=length(T1_intp);Tsamp=T1_intp(2)-T1_intp(1);fs=1/Tsamp;f1=(0:1/(N1-1):1)*fs;
startIdx = 3;
endIdx = round(0.04/(fs/N1))+1;
idx = startIdx:endIdx;
f_subset = f1(idx);
figure(20);clf;sgtitle('Patient 17');



%%
R_dft1 = R_dft1_patient17;
R_dft2 = R_dft2_patient17;
R_dft3 = R_dft3_patient17;
R_dft4 = R_dft4_patient17;
C_dft1 = C_dft1_patient17;
C_dft2 = C_dft2_patient17;
C_dft3 = C_dft3_patient17;
C_dft4 = C_dft4_patient17;

%% RESISTANCE
figure();
% Extract elements from 2 (without DC component) to endIdx
R1_dft_subset = abs(R_dft1(idx));
R2_dft_subset = abs(R_dft2(idx));
R3_dft_subset = abs(R_dft3(idx));
R4_dft_subset = abs(R_dft4(idx));

% find average
R1_dft_subset_mean_17 = mean(R1_dft_subset); 
R2_dft_subset_mean_17 = mean(R2_dft_subset); 
R3_dft_subset_mean_17 = mean(R3_dft_subset); 
R4_dft_subset_mean_17 = mean(R4_dft_subset);

R1_dft_subset_var_17 = mad(R1_dft_subset);
R2_dft_subset_var_17 = mad(R2_dft_subset);
R3_dft_subset_var_17 = mad(R3_dft_subset);
R4_dft_subset_var_17 = mad(R4_dft_subset);

R_mean = [R1_dft_subset_mean, R2_dft_subset_mean, R3_dft_subset_mean, R4_dft_subset_mean];
R_var = [R1_dft_subset_var, R2_dft_subset_var, R3_dft_subset_var, R4_dft_subset_var];

plot(f_subset,R1_dft_subset, '-+', f_subset,R2_dft_subset, '-o' , f_subset,R3_dft_subset, '-*', f_subset,R4_dft_subset,'-x');
legend('19.01.2019','20.01.2019','21.01.2019','23.01.2019'); xlabel('Frequency [Hz]'); ylabel('(1)');
%title('Relative resistance amplitude'); grid();

%%

figure();
errorbar((1:4),R_mean, R_var, '-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
%title('Mean amplitude (~20-240sec)');grid();
xlim([0.75,4.25]);xticklabels({'19.01','','20.01','','21.01','','23.01'}); xlabel('Date of measurement'); ylabel('(1)');

%% COMPLIANCE
% Extract elements from 2 (without DC component) to endIdx
figure();
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

plot(f_subset,C1_dft_subset, '-+', f_subset,C2_dft_subset, '-o' , f_subset,C3_dft_subset, '-*', f_subset,C4_dft_subset,'-x');
legend('19.01.2019','20.01.2019','21.01.2019','23.01.2019'); xlabel('Frequency [Hz]'); ylabel('(1)');
%title('Relative resistance amplitude'); grid();

%%

figure()
errorbar((1:4),C_mean, C_var, '-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
xlim([0.75,4.25]);xticklabels({'19.01','','20.01','','21.01','','23.01'}); xlabel('Date of measurement'); ylabel('(1)');
 