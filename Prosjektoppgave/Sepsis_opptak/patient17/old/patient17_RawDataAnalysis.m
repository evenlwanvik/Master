%% Analyse av pasient nr 17
% We have been given two snippets of measurements on different dates: 
% * 17/01/2019 (1) and 23/01/2019

load 20190117T145728_IQ_Sepsis-4min_traces;
Tmax1 = Tmax; Tmean1 = Tmean; Tmin1 = Tmin; Ts1 = Ts;
load 20190123T105641_IQ_Sepsis-4min_traces;
Tmax2 = Tmax; Tmean2 = Tmean; Tmin2 = Tmin; Ts2 = Ts;

%% Normal plot of both datasets

figure();clf;
subplot(2,2,1);plot(Ts1.t, Ts1.velocity);grid;
title("Raw velocity measurement - 20190117")
subplot(2,2,2);plot(Ts1.t, Ts1.ART);grid;
title("Raw pressure measurement - 20190117")
subplot(2,2,3);plot(Ts2.t, Ts2.velocity);grid;
title("Raw velocity measurement - 20190123")
subplot(2,2,4);plot(Ts2.t, Ts2.ART);grid;
title("Raw pressure measurement - 20190123")

%% Remove DC component and do frequency analysis


plot_dft(Ts1.velocity, Ts1.t, 10, '20190117 - Velocity DFT analysis', 1);
plot_dft(Ts1.ART, Ts1.t, 10, '20190117 - Pressure DFT analysis', 1);
plot_dft(Ts2.velocity, Ts2.t, 10, '20190123 - Velocity DFT analysis', 1);
plot_dft(Ts2.ART, Ts2.t, 10, '20190123 - Pressure DFT analysis', 1);

%%
figure();clf;

Fs = 100;
f1 = 10;f2 = 5;f3 = 15;
N = 2000
test_t = 0:1/Fs:(N-1)/Fs;
test_signal = 0.5*sin(2*pi*f1*test_t) + 2*sin(2*pi*f2*test_t) + sin(2*pi*f3*test_t);

%data = test_signal;
%fs = 1/(test_t(2)-test_t(1));
data = Ts2.velocity;
fs = 1/(Ts2.t(2)-Ts2.t(1));

N = length(data)

f = 0:Fs/(N-1):Fs;

windowed = hamming(N).*data;
DATA = mag2db( abs( fft(windowed) ) ) ;
semilogx(f, DATA);
grid();
title('Test signal')
ylabel('test signal dB')
xlabel('frequency')


%% Create a test signal and test our DFT function

% Test linearity

% DFT of the unit impulse
Fs = 100;
f1 = 10;
f2 = 5;
f3 = 15;
N = 2048
t = 0:1/Fs:(N-1)/Fs;
test = 0.5*sin(2*pi*f1*t) + 2*sin(2*pi*f2*t) + sin(2*pi*f3*t);
NFFT = N*10%2^nextpow2(N)
windowed = hamming(N)'.*test;
TEST = mag2db(abs(fft(windowed, NFFT)));
plot(TEST)


%% Test Windkessel model

%figure();clf;

% VELOCITY
T = Ts2.t(2)-Ts2.t(1)
fs = 1/T;
startIdx = round(8*fs)
endIdx = round(14*fs)
v = Ts2.velocity(startIdx:endIdx);
N = length(v)
subplot(4,1,1); plot(v); grid();

% PRESSURE
T = Ts2.t(2)-Ts2.t(1)
fs = 1/T;
startIdx = round(8*fs)
endIdx = round(14*fs)
p = Ts2.ART(startIdx:endIdx);
N = length(p)
subplot(4,1,2); plot(p); grid();

% WINDKESSEL model for finding velocity
P=fft(p);% fouriertransform av trykk p
w=2*pi*fs;%vinkelfrekvens
R=1; C=0.001;
Adm= 1/R - i*w*C;% admitans for paralellkobling av R og C
Z=1 ./Adm;% impedans
Q=Adm.*P;% flow i frekvensplan.  Fra Ohms lov P=Z*Q
q=real(ifft(Q));% flow i tidsplan
subplot(4,1,3); semilogx(abs(P)); grid();
subplot(4,1,4); plot(q); grid();



%% Functions

function plot_dft(data, t, N_segments, plot_title, semilog)
    % N_segments : Number of segments for average DFT 
    % data       : the dataset (made for raw data Ts)
    
    figure();clf;
    
    N = length(t); % Length of total measurement
    dt = t(2)-t(1);
    Fs = 1/dt
    
    % ### Remove offset
    DC = mean(data);
    data = data - DC;

    % ### Split dataset into multiple smaller datasets
    % ### e.g. binary scale 256 or 516 samples long
    % Set number of datapoints per segment (floored), i.e. we might be missing a few of the last datapoints
    segment_length = floor( N/N_segments ) ;       
    % Remove last ~0-100 samples of dataset and time such that we can divide it into equal chunks
    data_cut = data(1 : segment_length*N_segments);
    % Create dataset object (ds) and insert the data and time subarrays
    data_subarray = reshape(data_cut, segment_length, []);
    % find closest base-2 size for the DFT
    NFFT = 2^nextpow2(segment_length);
    % Perform a Hamming window, DFT the datasets and get it in decibels
    for i = 1:N_segments  
        windowed = hamming(segment_length).*data_subarray(:,i);
        plot(windowed)
        Vdb(:,i) = mag2db( abs( fft(windowed, NFFT) ) );
    end
    % Calculate the averaged DFT, the "mean()" function averaged in wrong
    % dimensions.. creating my own loop
    Vdb_mean = zeros(NFFT,1);
    for i = 1:NFFT 
        Vdb_mean(i) = mean(Vdb(i,:));
    end
    % Find the common frequency axis
    f = 0:Fs/(NFFT-1):Fs;
    % Plot fft for a random segment vs averaged fft
    if semilog == 1 
        semilogx(f, Vdb(:,2), f, Vdb_mean); ylabel('Power dB');
    else
        plot(f, Vdb(:,2), f, Vdb_mean); ylabel('Power');
    end
    title(plot_title); grid(); xlim([0 Fs]); xlabel('frequency'); 
end





