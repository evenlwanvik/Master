%% Extract parameters over time for two-element Windkessel model.

addpath('Master/Prosjektoppgave/Sepsis_opptak/patient17/data') 
%load 20190117T145728_IQ_Sepsis-4min_traces;
%date = '17.01.2019';
load 20190123T105641_IQ_Sepsis-4min_traces;
date = '23.01.2019';
name = append('Sepsis 4min_traces ', date)

%% Analyse measurements and connect end-points to avoid FFT leakage

delay=0.07;
tED=Tmean.tED-delay;
Ncycles=length(tED);

% Get data
n=140; % what cycle to extract
tIdx = find( Ts.t>tED(n) & Ts.t<tED(n+20)); % all samples within tED1 - tED2 window
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
v_after = v - v_line'+v(1); %v_after = v_after-mean(v_after);
p_after = p - p_line'+p(1); %p_after = p_after-mean(p_after);

% Plot after
subplot(2,2,3);plot(t,v_after,x_range,[v_after(1),v_after(end)]);title('Velocity after');xlim(x_range);grid();
subplot(2,2,4);plot(t,p_after,x_range,[p_after(1),p_after(end)]);title('Pressure after');xlim(x_range);grid();

%% Analyse the frequency spectra of the measurements before and after

f = (0:1/(N-1):1)*fs;
%V = fft(hamming(N,'periodic').*v); V_after = fft(hamming(N,'periodic').*v_after);
%P = fft(hamming(N,'periodic').*p); P_after = fft(hamming(N,'periodic').*p_after);
V = fft(v); V_after = fft(v_after);
P = fft(p); P_after = fft(p_after);

% plot all
figure();clf;sgtitle(append('V and P DFT - ', date));
subplot(2,2,1);plot(f,abs(V));title('Velocity before');grid();xlim([0 500]);ylim([0, 0.1]);
subplot(2,2,2);plot(f,abs(P));title('Pressure before');grid();xlim([0 500]);ylim([0, 10]);
subplot(2,2,3);plot(f,abs(V_after));title('Velocity after');grid();xlim([0 500]);ylim([0, 0.1]);
subplot(2,2,4);plot(f,abs(P_after));title('Pressure after');grid();xlim([0 500]);ylim([0, 10]);

%% Analyse the frequency spectra of the impedance before and after

w=2*pi*f;

% Get the index and 1. harmonic frequency
[Av, i_Fv] = max(abs(V(2:N-1))); i_Fv=i_Fv+1;

% Before
Z=P./V; 
R=Z(1);
C = abs( 1i*(1/Z(i_Fv)-1/R) / (2*pi*f(i_Fv)) )
%C = 1e-3;
Adm= 1/R - 1i*w*C;
Z_WK=1./Adm;

% After
Z_after=P_after./V_after;
R_after=Z_after(1);
C_after = abs( 1i*(1/Z_after(i_Fv)-1/R_after) / (2*pi*f(i_Fv)) );
%C_after = 1e-3;
Adm_WK_after= 1/R_after - 1i*w*C_after;
Z_WK_after=1./Adm_WK_after;

% plot
figure();clf;sgtitle(append('Z DFT - ', date));
subplot(2,1,1);plot(f,abs(Z),f,abs(Z_WK));title('Z before');grid();xlim([0,500]);ylim([0,2e4])
subplot(2,1,2);plot(f,abs(Z_after),f,abs(Z_WK_after));title('Z after');grid();ylim([0,2e4])



%% Create plot of compliance and resistance

[R_arr, C_arr, t_axis] = extractParameters(Ts, Tmean);

% Average?
%R_arr = movmean(R_arr, 3); C_arr = movmean(C_arr, 3);

figure();clf;title(append('Normal - Resistance (l) and compliance (r) ', date));
yyaxis left
plot(R_arr); xlim([0, 316]); ylim([0, 7000])
yyaxis right
plot(C_arr); ylim([0, 1.2e-3]);grid();


%% Functions



% The two following functions extracts the parameters
function [R_arr, C_arr, t_axis] = extractParameters(Ts, Tmean)
    delay = 0.07;
    tED = Tmean.tED-delay;
    nBeats = 20; % number of beats we analyse per "cycle"
    totCycles = floor(length(tED/nBeats));
    Tsamp = tED(2)-tED(1); Fsamp = 1/Tsamp;
    C_arr = []; R_arr = []; t_axis = [];
    for k = 1:totCycles-nBeats
        tIdx = find( Ts.t>tED(k) & Ts.t<tED(k+nBeats)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx); N=length(t);
        % Create line between start and end of signal
        v_line = linspace(v(1), v(end), N);
        p_line = linspace(p(1), p(end), N);
        v_after = v - v_line'+v(1); 
        p_after = p - p_line'+p(1); 
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
    
    [Av, i_Fv] = max(abs(V(2:N-1))); i_Fv=i_Fv+1;
    C = abs( 1i*(1/Z(i_Fv)-1/R) / (2*pi*f(i_Fv)) );
    %C = abs( ( 1/R - 1/Z(2) ) / (2*pi*f(2)) );

end