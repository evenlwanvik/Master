%% Extract parameters over time for two-element Windkessel model.

%load 20190117T145728_IQ_Sepsis-4min_traces;
%date = '17.01.2019'
load 20190123T105641_IQ_Sepsis-4min_traces;
date = '23.01.2019'

name = append('Sepsis 4min_traces ', date)

figure();clf;sgtitle(name);
[R_arr, C_arr, z_WK, t_axis] = loop_through_dataset(Ts, Tmean);

% Average?
R_arr = movmean(R_arr, 3); C_arr = movmean(C_arr, 3);

figure();clf;title(append('Resistance (l) and compliance (r) ', date));
yyaxis left
plot(t_axis, R_arr)
ylim([0, 1500])
yyaxis right
plot(t_axis, C_arr)


%% Reconstruct impedance with two-element Windkessel model
figure();clf;sgtitle(name);

delay = 0.07;
tED = Tmean.tED-delay;
    
Z_WK = abs(fft(z_WK));
T = Z_WK();
plot(abs(Z_WK));


%% functions 

function [R_arr, C_arr, Z_WK_arr, t_axis] = loop_through_dataset(Ts, Tmean)
    delay = 0.07;
    tED = Tmean.tED-delay;
    totCycles = length(tED)
    Tsamp = tED(2)-tED(1); Fsamp = 1/Tsamp
    C_arr = []; R_arr = []; Z_WK_arr = []; t_axis = [];

    % Create axes for continious plotting
    ax1 = subplot(2,2,1); ax2 = subplot(2,2,2); ax3 = subplot(2,2,3); ax4 = subplot(2,2,4);
    hold(ax1, 'on'); hold(ax2, 'on'); hold(ax3, 'on'); hold(ax4, 'on');
    for k = 1:totCycles-1
        tIdx = find( Ts.t>tED(k) & Ts.t<tED(k+1)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); 
        v=Ts.velocity(tIdx);
        p=Ts.ART(tIdx);       
        [C, R, Z_WK] = get_parameters(t, p, v, ax1, ax2, ax3, ax4);
        R_arr = [R_arr, R];
        C_arr = [C_arr, C];
        Z_WK_arr = [Z_WK_arr, Z_WK];
        t_axis = [t_axis, tED(k)];
    end
    %hold off;
end

function [C, R, Z_WK] = get_parameters(t, p, v, ax1, ax2, ax3, ax4)

    Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);
    f = (0:1/(N-1):1)*fs;
    V = fft(v); 
    P = fft(p); 
    Z = P./V;
    R = Z(1); % zero freq component is the real value of Z

    % WINDKESSEL model for finding velocity
    % I know that second freq index is the first harmonic of the heart cycle
    C = abs( ( 1/R - 1/Z(2) ) / (2*pi*f(2)) );
    Adm = 1/R + i*2*pi*f(2)*C;
    Q = Adm.*P;
    q=abs(real(ifft(Q)));% flow i tidsplan
    Z_WK = 1/Adm;
    % Continious plot of inverse
    %axes(ax1); plot(t,p);title('Heart cycle - pressure')
    %axes(ax2); plot(t,v);title('Heart cycle - velocity')
    %axes(ax3); plot(t,q);title('Velocity remade by IFFT')
    %axes(ax4); plot(t,abs(Z_WK)); title('Windkessel impedance')
end

