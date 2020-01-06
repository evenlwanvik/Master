% Run "week_44/compliance_patient##.m first"

% we wants oscillations from 20 sec to 2min, i.e. 0-0.05  Hz
N1=length(T1_intp);Tsamp=T1_intp(2)-T1_intp(1);fs=1/Tsamp;f1=(0:1/(N1-1):1)*fs;
startIdx = 3;
endIdx = round(0.04/(fs/N1))+1;
idx = startIdx:endIdx;
f_subset = f1(idx);

%%
noradrenaline30 = [.6, .32, .1];
noradrenaline29 = [.34, .3, .3, .26, .26, .16, .04];
noradrenaline28 = [.2, .16, 0.0];
noradrenaline24 = [.28, .1, .02, .02];
noradrenaline19 = [.11, .11, .1, .6];
noradrenaline18 = [0.07, 0.12, 0.13, 0.11];
noradrenaline17 = [.17, .12, 0.01, 0];
noradrenaline16 = [.33, .3, .33, .04];
noradrenaline15 = [.16, .02, .18, .24];

%%

figure(15);

% Patient 15
R_mean_15 = [R1_dft_subset_mean_15, R2_dft_subset_mean_15, R3_dft_subset_mean_15, R4_dft_subset_mean_15];
R_var_15 = [R1_dft_subset_var_15, R2_dft_subset_var_15, R3_dft_subset_var_15, R4_dft_subset_var_15];
subplot(4,2,1); title('patient 15');
yyaxis left
errorbar((1:4),R_mean_15, R_var_15, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline15, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','21.01','23.01'}); 
ylim([0, 0.3]); 

% Patient 16
R_mean_16 = [R1_dft_subset_mean_16, R2_dft_subset_mean_16, R3_dft_subset_mean_16, R4_dft_subset_mean_16];
R_var_16 = [R1_dft_subset_var_16, R2_dft_subset_var_16, R3_dft_subset_var_16, R4_dft_subset_var_16];
subplot(4,2,2); title('patient 16');
yyaxis left
errorbar((1:4),R_mean_16, R_var_16, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline16, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','21.01','23.01'}); 
ylim([0, 0.40]); 

% Patient 18
R_mean_18 = [R1_dft_subset_mean_18, R2_dft_subset_mean_18, R3_dft_subset_mean_18, R4_dft_subset_mean_18];
R_var_18 = [R1_dft_subset_var_18, R2_dft_subset_var_18, R3_dft_subset_var_18, R4_dft_subset_var_18];
subplot(4,2,3); title('patient 18');
yyaxis left
errorbar((1:4),R_mean_18, R_var_18, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline18, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.15]); 

% Patient 19
R_mean_19 = [R1_dft_subset_mean_19, R2_dft_subset_mean_19, R3_dft_subset_mean_19, R4_dft_subset_mean_19];
R_var_19 = [R1_dft_subset_var_19, R2_dft_subset_var_19, R3_dft_subset_var_19, R4_dft_subset_var_19];
subplot(4,2,4); title('patient 19');
yyaxis left
errorbar((1:4),R_mean_19, R_var_19, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline19, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.65]); 

% Patient 24
R_mean_24 = [R1_dft_subset_mean_24, R2_dft_subset_mean_24, R3_dft_subset_mean_24, R4_dft_subset_mean_24];
R_var_24 = [R1_dft_subset_var_24, R2_dft_subset_var_24, R3_dft_subset_var_24, R4_dft_subset_var_24];
subplot(4,2,5); title('patient 24');
yyaxis left
errorbar((1:4),R_mean_24, R_var_24, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline24, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.35]); 

% Patient 28
R_mean_28 = [R1_dft_subset_mean_28, R2_dft_subset_mean_28, R3_dft_subset_mean_28];
R_var_28 = [R1_dft_subset_var_28, R2_dft_subset_var_28, R3_dft_subset_var_28];
subplot(4,2,6); title('patient 28');
yyaxis left
errorbar((1:3),R_mean_28, R_var_28, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:3), noradrenaline28, '--o', 'MarkerSize', 8); 
xlim([0.75,3.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.25]); 

% Patient 29
R_mean_29 = [R1_dft_subset_mean_29, R2_dft_subset_mean_29, R3_dft_subset_mean_29, R4_dft_subset_mean_29, R5_dft_subset_mean_29, R6_dft_subset_mean_29, R7_dft_subset_mean_29];
R_var_29 = [R1_dft_subset_var_29, R2_dft_subset_var_29, R3_dft_subset_var_29, R4_dft_subset_var_29, R5_dft_subset_var_29, R6_dft_subset_var_29, R7_dft_subset_var_29];
subplot(4,2,7); title('patient 29');
yyaxis left
errorbar((1:7), R_mean_29, R_var_29, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
ylim([0, 0.035]); 
yyaxis right
plot((1:7), noradrenaline29, '--o', 'MarkerSize', 8); 
xlim([0.75,7.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.35]); 

%% COMPLIANCE

figure(16);

% Patient 15
C_mean_15 = [C1_dft_subset_mean_15, C2_dft_subset_mean_15, C3_dft_subset_mean_15, C4_dft_subset_mean_15];
C_var_15 = [C1_dft_subset_var_15, C2_dft_subset_var_15, C3_dft_subset_var_15, C4_dft_subset_var_15];
subplot(4,2,1); title('patient 15');
yyaxis left
errorbar((1:4),C_mean_15, C_var_15, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline15, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','21.01','23.01'}); 
ylim([0, 0.3]); 

% Patient 16
C_mean_16 = [C1_dft_subset_mean_16, C2_dft_subset_mean_16, C3_dft_subset_mean_16, C4_dft_subset_mean_16];
C_var_16 = [C1_dft_subset_var_16, C2_dft_subset_var_16, C3_dft_subset_var_16, C4_dft_subset_var_16];
subplot(4,2,2); title('patient 16');
yyaxis left
errorbar((1:4),C_mean_16, C_var_16, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline16, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','21.01','23.01'}); 
ylim([0, 0.40]); 

% Patient 18
C_mean_18 = [C1_dft_subset_mean_18, C2_dft_subset_mean_18, C3_dft_subset_mean_18, C4_dft_subset_mean_18];
C_var_18 = [C1_dft_subset_var_18, C2_dft_subset_var_18, C3_dft_subset_var_18, C4_dft_subset_var_18];
subplot(4,2,3); title('patient 18');
yyaxis left
errorbar((1:4),C_mean_18, C_var_18, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline18, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.15]); 

% Patient 19
C_mean_19 = [C1_dft_subset_mean_19, C2_dft_subset_mean_19, C3_dft_subset_mean_19, C4_dft_subset_mean_19];
C_var_19 = [C1_dft_subset_var_19, C2_dft_subset_var_19, C3_dft_subset_var_19, C4_dft_subset_var_19];
subplot(4,2,4); title('patient 19');
yyaxis left
errorbar((1:4),C_mean_19, C_var_19, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline19, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.65]); 

% Patient 24
C_mean_24 = [C1_dft_subset_mean_24, C2_dft_subset_mean_24, C3_dft_subset_mean_24, C4_dft_subset_mean_24];
C_var_24 = [C1_dft_subset_var_24, C2_dft_subset_var_24, C3_dft_subset_var_24, C4_dft_subset_var_24];
subplot(4,2,5); title('patient 24');
yyaxis left
errorbar((1:4),C_mean_24, C_var_24, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:4), noradrenaline24, '--o', 'MarkerSize', 8); 
xlim([0.75,4.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.35]); 

% Patient 28
C_mean_28 = [C1_dft_subset_mean_28, C2_dft_subset_mean_28, C3_dft_subset_mean_28];
C_var_28 = [C1_dft_subset_var_28, C2_dft_subset_var_28, C3_dft_subset_var_28];
subplot(4,2,6); title('patient 28');
yyaxis left
errorbar((1:3),C_mean_28, C_var_28, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
yyaxis right
plot((1:3), noradrenaline28, '--o', 'MarkerSize', 8); 
xlim([0.75,3.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.25]); 

% Patient 29
C_mean_29 = [C1_dft_subset_mean_29, C2_dft_subset_mean_29, C3_dft_subset_mean_29, C4_dft_subset_mean_29, C5_dft_subset_mean_29, C6_dft_subset_mean_29, C7_dft_subset_mean_29];
C_var_29 = [C1_dft_subset_var_29, C2_dft_subset_var_29, C3_dft_subset_var_29, C4_dft_subset_var_29, C5_dft_subset_var_29, C6_dft_subset_var_29, C7_dft_subset_var_29];
subplot(4,2,7); title('patient 29');
yyaxis left
errorbar((1:7), C_mean_29, C_var_29, '-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',25)
ylim([0, 0.035]); 
yyaxis right
plot((1:7), noradrenaline29, '--o', 'MarkerSize', 8); 
xlim([0.75,7.25]);xticklabels({'19.01','20.01','23.01'}); 
ylim([0, 0.35]); 