
%%
% Extract elements from 2 (without DC component) to endIdx
C1_dft_subset = abs(C_dft1(3:7));
C2_dft_subset = abs(C_dft2(3:7));
C3_dft_subset = abs(C_dft3(3:7));
C4_dft_subset = abs(C_dft4(3:7));



% find average
C1_dft_subset_mean = mean(C_dft1(3:7));
C2_dft_subset_mean = mean(C_dft2(3:7));
C3_dft_subset_mean = mean(C_dft3(3:7));
C4_dft_subset_mean = mean(C_dft4(3:7));

C_mean_trend_patient18 = [C1_dft_subset_mean, C2_dft_subset_mean, C3_dft_subset_mean, C4_dft_subset_mean];


%%
R = {R_dft1_patient18, R_dft1_patient17, R_dft1_patient16, R_dft1_patient15;
    R_dft2_patient18, R_dft2_patient17, R_dft2_patient16, R_dft1_patient15;
    R_dft3_patient18, R_dft3_patient17, R_dft3_patient16, R_dft1_patient15;
    R_dft4_patient18, R_dft4_patient17, R_dft4_patient16, R_dft1_patient15;};
C = {C_dft1_patient18, C_dft1_patient17, C_dft1_patient16, C_dft1_patient15;
    C_dft2_patient18, C_dft2_patient17, C_dft2_patient16, C_dft1_patient15;
    C_dft3_patient18, C_dft3_patient17, C_dft3_patient16, C_dft1_patient15;
    C_dft4_patient18, C_dft4_patient17, C_dft4_patient16, C_dft1_patient15;};
figure()
subplot(2,2,1); xlim([0.75,4.25]); xticks([1,2,3,4]); xticklabels({'11.12', '12.12', '18.12', '19.12'})
yyaxis left
plot_mean_freq(R_dft1_patient18, R_dft2_patient18, R_dft3_patient18, R_dft4_patient18);ylabel('$\hat{R}(f)$','Interpreter','latex')
yyaxis right
plot_mean_freq(C_dft1_patient18, C_dft2_patient18, C_dft3_patient18, C_dft4_patient18);ylabel('$\hat{C}(f)$','Interpreter','latex')
subplot(2,2,2); xlim([0.75,4.25]); xticks([1,2,3,4]); xticklabels({'13.12', '13.12', '14.12', '17.12'})
yyaxis left
plot_mean_freq(R_dft1_patient17, R_dft2_patient17, R_dft3_patient17, R_dft4_patient17);ylabel('$\hat{R}(f)$','Interpreter','latex')
yyaxis right
plot_mean_freq(C_dft1_patient17, C_dft2_patient17, C_dft3_patient17, C_dft4_patient17);ylabel('$\hat{C}(f)$','Interpreter','latex')
subplot(2,2,3); xlim([0.75,4.25]); xticks([1,2,3,4]); xticklabels({'19.01', '20.01', '21.01', '23.01'})
yyaxis left
plot_mean_freq(R_dft1_patient16, R_dft2_patient16, R_dft3_patient16, R_dft4_patient16);ylabel('$\hat{R}(f)$','Interpreter','latex')
yyaxis right
plot_mean_freq(C_dft1_patient16, C_dft2_patient16, C_dft3_patient16, C_dft4_patient16);ylabel('$\hat{C}(f)$','Interpreter','latex')
subplot(2,2,4); xlim([0.75,4.25]); xticks([1,2,3,4]); xticklabels({'19.01', '20.01', '21.01', '23.01'})
yyaxis left
plot_mean_freq(R_dft1_patient15, R_dft2_patient15, R_dft3_patient15, R_dft4_patient15);ylabel('$\hat{R}(f)$','Interpreter','latex')
yyaxis right
plot_mean_freq(C_dft1_patient15, C_dft2_patient15, C_dft3_patient15, C_dft4_patient15);ylabel('$\hat{C}(f)$','Interpreter','latex')
%xticks([1,2,3,4]); xticklabels({'first', 'second', 'third', 'fourth'})

%%



figure()
hold on
plot_mean_freq(R_dft1_patient18, R_dft2_patient18, R_dft3_patient18, R_dft4_patient18);
plot_mean_freq(R_dft1_patient17, R_dft2_patient17, R_dft3_patient17, R_dft4_patient17);
plot_mean_freq(R_dft1_patient16, R_dft2_patient16, R_dft3_patient16, R_dft4_patient16);
plot_mean_freq(R_dft1_patient15, R_dft2_patient15, R_dft3_patient15, R_dft4_patient15);
hold off
xlim([0.75,4.25]);xticks([1,2,3,4]);
xticklabels({'first', 'second', 'third', 'fourth'})
legend('Patient 18', 'Patient 17', 'Patient 16', 'Patient 15')

%%

function [] = plot_mean_freq(dft1, dft2, dft3, dft4)
    mean1 = mean(abs(dft1(3:7)));
    mean2 = mean(abs(dft2(3:7)));
    mean3 = mean(abs(dft3(3:7)));
    mean4 = mean(abs(dft4(3:7)));
    trend = [mean1, mean2, mean3, mean4];
    plot(1:4, trend)
end