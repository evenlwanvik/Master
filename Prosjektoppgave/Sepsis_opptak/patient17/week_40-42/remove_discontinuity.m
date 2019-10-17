%% Extract parameters over time for two-element Windkessel model.

addpath('Master/Prosjektoppgave/Sepsis_opptak/patient17/data') 
load 20190117T145728_IQ_Sepsis-4min_traces;
date = '17.01.2019';
%load 20190123T105641_IQ_Sepsis-4min_traces;
%date = '23.01.2019';
name = append('Sepsis 4min_traces ', date)

%% Analyse measurements and connect end-points to avoid FFT leakage

delay=0.07;
tED=Tmean.tED-delay;
Ncycles=length(tED);

% Get data
n=140; % what cycle to extract
tIdx = find( Ts.t>tED(n) & Ts.t<tED(n+1)); % all samples within tED1 - tED2 window
t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx);
Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);

% Plot before
figure();clf;sgtitle(append('Removing discontinuities - measurements - ', date));
x_range = [t(1),t(end)];
subplot(2,2,1);plot(t,v,x_range,[v(1),v(end)]);title('Velocity before');xlim(x_range);grid();
subplot(2,2,2);plot(t,p,x_range,[p(1),p(end)]);title('Pressure before');xlim(x_range);grid();

% Create line between start and end of signal
v_line = linspace(v(1), v(end), N);
p_line = linspace(p(1), p(end), N);
v_after = v - v_line'+v(1);
p_after = p - p_line'+p(1);

% Plot after
subplot(2,2,3);plot(t,v_after,x_range,[v_after(1),v_after(end)]);title('Velocity after');xlim(x_range);grid();
subplot(2,2,4);plot(t,p_after,x_range,[p_after(1),p_after(end)]);title('Pressure after');xlim(x_range);grid();

%% Filter signal?
% Filter signal
[zeros,poles,gain] = butter(5,10/fs,'low');
sos = zp2sos(zeros,poles,gain); % sos stands for second-order biquadric IIR digital filtering
v_filt = sosfilt(sos, v_after); v_filt = v_filt - linspace(v_filt(1), v_filt(end), N)';
p_filt = sosfilt(sos, p_after); p_filt = p_filt - linspace(p_filt(1), p_filt(end), N)';

% Plot filtered
subplot(3,2,5);plot(t,v_filt,x_range,[v_filt(1),v_filt(end)]);title('Velocity filtered');xlim(x_range);grid();
subplot(3,2,6);plot(t,p_filt,x_range,[p_filt(1),p_filt(end)]);title('Pressure filtered');xlim(x_range);grid();
%% Analyse the frequency spectra of the measurements before and after

f = (0:1/(N-1):1)*fs;
V = fft(v); V_after = fft(v_after);
P = fft(p); P_after = fft(p_after);

% plot all
figure();clf;sgtitle(append('V and P DFT - ', date));
subplot(2,2,1);plot(f,mag2db(abs(V)));title('Velocity before');grid();%xlim([0 500]);ylim([0, 0.2]);
subplot(2,2,2);plot(f,mag2db(abs(P)));title('Pressure before');grid();%xlim([0 500]);ylim([0, 200]);
subplot(2,2,3);plot(f,mag2db(abs(V_after)));title('Velocity after');grid();%xlim([0 500]);ylim([0, 0.2]);
subplot(2,2,4);plot(f,mag2db(abs(P_after)));title('Pressure after');grid();%xlim([0 500]);ylim([0, 200]);

%% Analyse the frequency spectra of the impedance before and after

w=2*pi*f;

% Before
Z=P./V; 
R=Z(1);
C = abs( 1i*(1/Z(2)-1/R) / (2*pi*f(2)) );
Adm= 1/R - 1i*w*C;
Z_WK=1./Adm;

% After
Z_after=P_after./V_after;
R_after=Z_after(1);
C_after = abs( 1i*(1/Z_after(2)-1/R_after) / (2*pi*f(2)) );
Adm_WK_after= 1/R_after - 1i*w*C_after;
Z_WK_after=1./Adm_WK_after;

% plot
figure();clf;sgtitle(append('Z DFT - ', date));
subplot(2,1,1);plot(f,abs(Z),f,abs(Z_WK));title('Z before');grid();xlim([0,100])
subplot(2,1,2);plot(f,abs(Z_after),f,abs(Z_WK_after));title('Z after');grid();xlim([0,100])

%% Do a average of all 

[Z_before_arr, Z_after_arr] = average_impedance_spectra(Ts, Tmean);
Z_before_mean = mean(Z_before_arr, 2);
Z_after_mean = mean(Z_after_arr, 2);

Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(Z_before_mean);
f = (0:1/(N-1):1)*fs;

% plot
figure();clf;sgtitle(append('Z DFT averaged - ', date));
subplot(2,1,1);plot(f,abs(Z_before_mean));title('Z before');grid();xlim([0,500]);ylim([0,7000]);
subplot(2,1,2);plot(f,abs(Z_after_mean));title('Z after');grid();xlim([0,500]);ylim([0,7000]);

%% Test resonant peaks in averaged impedance spectra

i=10;x=0;
resonant_peaks=[];
%window=[140, 150];
%fIdx = find(f>window(1) & f<window(2))
%[val, idx] = max(Z_before_mean(fIdx));
%resonant_peak = f(idx+fIdx(1)-1)
resonant_peak = [];
while x<150
    x = f(i);
    if Z_before_mean(i)>4000 resonant_peaks= [resonant_peaks, f(i)]; end
    i = i+1;
end
resonant_peaks

%% Create plot of compliance and resistance

[R_arr, C_arr, t_axis] = extractParameters(Ts, Tmean);

% Average?
%R_arr = movmean(R_arr, 3); C_arr = movmean(C_arr, 3);

figure();clf;title(append('Resistance (l) and compliance (r) ', date));
yyaxis left
plot(R_arr); xlim([0, 316]); ylim([0, 9000])
yyaxis right
plot(C_arr); ylim([0, 2.6e-4]); 

%% analyse all cycles to check what our main characteristics are

function all_cycles()
    figure();clf;sgtitle(append('Heart cycle analysis - ', date));
    for i=1:Ncycles
        % Get data
        tIdx = find( Ts.t>tED(i) & Ts.t<tED(i+1)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx);

        % Plot before

        subplot(1,2,1);plot(t,v,[t(1),t(end)],[v(1),v(end)]);title('Heart cycle - velocity');grid();
        subplot(1,2,2);plot(t,p,[t(1),t(end)],[p(1),p(end)]);title('Heart cycle - pressure');grid();

        % Plot after

        pause(1);
    end
end

function analyze_impedance(Ts, Tmean)
    delay=0.07;
    tED=Tmean.tED-delay;
    Ncycles=length(tED);
    figure();clf;sgtitle(append('Z DFT - ', date));
    for i=1:Ncycles
        tIdx = find( Ts.t>tED(i) & Ts.t<tED(i+1)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx);
        Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);
        % Create line between start and end of signal
        v_line = linspace(v(1), v(end), N);
        p_line = linspace(p(1), p(end), N);
        v_after = v - v_line'; %v_after = v_after-mean(v_after);
        p_after = p - p_line'; %p_after = p_after-mean(p_after);
        % V and P freq spectra
        f = (0:1/(N-1):1)*fs; w=2*pi*f;
        V = fft(v); V_after = fft(v_after);
        P = fft(p); P_after = fft(p_after);

        % Z Before
        Z=P./V; R=Z(1); C = abs( 1i*(1/Z(2)-1/R) / (2*pi*f(2)) );
        Adm= 1/R - 1i*w*C; Z_WK=1./Adm;

        % Z After
        Z_after=P_after./V_after; R_after=Z_after(1);
        C_after = abs( 1i*(1/Z_after(2)-1/R_after) / (2*pi*f(2)) );
        Adm_WK_after= 1/R_after - 1i*w*C_after; Z_WK_after=1./Adm_WK_after;

        % plot
        subplot(2,1,1);plot(f,abs(Z),f,abs(Z_WK));title('Z before');grid();xlim([0,100]);ylim([0,3000]);
        subplot(2,1,2);plot(f,abs(Z_after),f,abs(Z_WK_after));title('Z after');grid();xlim([0,100]);ylim([0,3000]);
        
        pause(1);
    end    
end

function [Z_before_arr, Z_after_arr] = average_impedance_spectra(Ts, Tmean)
    delay=0.05;
    tED=Tmean.tED-delay;
    Ncycles=length(tED);
	% Just use a random size of a cycle for averaging dimensions
    %Z_before_arr = zeros(length(find( Ts.t>tED(15) & Ts.t<tED(15+1))), Ncycles);
    %Z_after_arr = Z_before;
    for i=1:Ncycles
        if i>=length(tED) continue; end
        tIdx = find( Ts.t>tED(i) & Ts.t<tED(i+1)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx);
        Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);
        % Create line between start and end of signal
        v_line = linspace(v(1), v(end), N);
        p_line = linspace(p(1), p(end), N);
        v_after = v - v_line'; %v_after = v_after-mean(v_after);
        p_after = p - p_line'; %p_after = p_after-mean(p_after);
        % V and P freq spectra
        f = (0:1/(N-1):1)*fs; w=2*pi*f;
        V = fft(v); V_after = fft(v_after);
        P = fft(p); P_after = fft(p_after);

        % Z Before
        Z=P./V; R=Z(1); C = abs( 1i*(1/Z(2)-1/R) / (2*pi*f(2)) );
        Adm= 1/R - 1i*w*C; Z_WK=1./Adm;

        % Z After
        Z_after=P_after./V_after; R_after=Z_after(1);
        C_after = abs( 1i*(1/Z_after(2)-1/R_after) / (2*pi*f(2)) );
        Adm_WK_after= 1/R_after - 1i*w*C_after; Z_WK_after=1./Adm_WK_after;

        % plot
        for k=1:313 % found that all of the cycles were at least 314 indeces long
            Z_before_arr(k,i) = abs(Z(k));
            Z_after_arr(k,i) = abs(Z_after(k));
        end

    end
end

function [R_arr, C_arr, t_axis] = extractParameters(Ts, Tmean)
    delay = 0.07;
    tED = Tmean.tED-delay;
    totCycles = length(tED);
    Tsamp = tED(2)-tED(1); Fsamp = 1/Tsamp;
    C_arr = []; R_arr = []; t_axis = [];
    for k = 1:totCycles-1
        tIdx = find( Ts.t>tED(k) & Ts.t<tED(k+1)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx); N=length(t);
        % Create line between start and end of signal
        v_line = linspace(v(1), v(end), N);
        p_line = linspace(p(1), p(end), N);
        v_after = v - v_line'; %v_after = v_after-mean(v_after);
        p_after = p - p_line'; %p_after = p_after-mean(p_after);
        % get parameters
        [C, R] = get_parameters(t, p_after, v_after);
        R_arr = [R_arr, R]; C_arr = [C_arr, C]; t_axis = [t_axis, tED(k)];
    end
end

function [C, R] = get_parameters(t, p, v)

    Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);
    f = (0:1/(N-1):1)*fs;
    V = fft(v); P = fft(p); Z = P./V;
    R = Z(1); % zero freq component is the real value of Z
    C = abs( ( 1/R - 1/Z(2) ) / (2*pi*f(2)) );

end
