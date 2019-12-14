
%figure(12);clf;sgtitle('Patient 18 - Compliance');
subplot(4,1,1);
plot(T1_intp,C1_intp);ylim([3e-4,8e-4]);ylabel('Compliance');title(name1);
subplot(4,1,2);
plot(T2_intp,C2_intp);ylim([4e-4,8e-4]);ylabel('Compliance');title(name2);
subplot(4,1,3);
plot(T3_intp,C3_intp);ylim([5e-4,8e-4]);ylabel('Compliance');title(name3);
subplot(4,1,4);xlabel('Time [s]');
plot(T4_intp,C4_intp);ylim([4.5e-4,7.5e-4]);ylabel('Compliance');title(name4);
xlabel('Time [s]')

%%

subplot(4,1,1);
plot(T1_intp,R1_intp);ylim([3000,6000]);ylabel('Resistance');title(name1);
subplot(4,1,2);
plot(T2_intp,R2_intp);ylim([1500,3500]);ylabel('Resistance'); title(name2);
subplot(4,1,3);
plot(T3_intp,R3_intp);ylim([1000,3000]);ylabel('Resistance'); title(name3);
subplot(4,1,4);xlabel('Time [s]');
plot(T4_intp,R4_intp);ylim([900,1400]);ylabel('Resistance');title(name4);
