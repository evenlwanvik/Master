
figure(4)
yyaxis left
plot(T1_intp,R1_intp);ylim([900,1400]);ylabel('Resistance')
yyaxis right
plot(T1_intp,C1_intp);ylim([1e-4,9e-4]);ylabel('Compliance')
title(name4);xlabel('Time [s]');


%%

subplot(4,1,1);
plot(T1_intp,R1_intp);ylim([3000,6000]);ylabel('Resistance');title(name1);
subplot(4,1,2);
plot(T2_intp,R2_intp);ylim([1500,3500]);ylabel('Resistance'); title(name2);
subplot(4,1,3);
plot(T3_intp,R3_intp);ylim([1000,3000]);ylabel('Resistance'); title(name3);
subplot(4,1,4);xlabel('Time [s]');
plot(T4_intp,R4_intp);ylim([900,1400]);ylabel('Resistance');title(name4);
xlabel('Time [s]')

%%

C_dft1_patient18 = C_dft1;
C_dft2_patient18 = C_dft2;
C_dft3_patient18 = C_dft3;
C_dft4_patient18 = C_dft4;
R_dft1_patient18 = R_dft1;
R_dft2_patient18 = R_dft2;
R_dft3_patient18 = R_dft3;
R_dft4_patient18 = R_dft4;

%%
C_dft1_patient17 = C_dft1;
C_dft2_patient17 = C_dft2;
C_dft3_patient17 = C_dft3;
C_dft4_patient17 = C_dft4;
R_dft1_patient17 = R_dft1;
R_dft2_patient17 = R_dft2;
R_dft3_patient17 = R_dft3;
R_dft4_patient17 = R_dft4;

%%
C_dft1_patient16 = C_dft1;
C_dft2_patient16 = C_dft2;
C_dft3_patient16 = C_dft3;
C_dft4_patient16 = C_dft4;
R_dft1_patient16 = R_dft1;
R_dft2_patient16 = R_dft2;
R_dft3_patient16 = R_dft3;
R_dft4_patient16 = R_dft4;

%%
C_dft1_patient15 = C_dft1;
C_dft2_patient15 = C_dft2;
C_dft3_patient15 = C_dft3;
C_dft4_patient15 = C_dft4;
R_dft1_patient15 = R_dft1;
R_dft2_patient15 = R_dft2;
R_dft3_patient15 = R_dft3;
R_dft4_patient15 = R_dft4;