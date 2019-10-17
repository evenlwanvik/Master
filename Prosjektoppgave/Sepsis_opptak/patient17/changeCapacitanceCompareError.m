%% Compare the first 5 Hz of actual impedance and WK model impedance

load 20190117T145728_IQ_Sepsis-4min_traces;
date = '17.01.2019'
%load 20190123T105641_IQ_Sepsis-4min_traces;
%date = '23.01.2019'
name = append('Sepsis 4min_traces ', date)

err = []; err4Freq = []; 
C = 0.00001:0.00001:0.001;
cycles= 1:10:310; % we got about 330-350 timestamps for each cycle
for C_k=C
    errSum = 0;
    err4FreqSum = 0;
    for cycleStart=cycles
        
        %-------- measured data from single heart cycle --------
        delay = 0.07;
        tED = Tmean.tED-delay;
        %cycleStart=33; % choose heart cycle (~300 possibilities)
        nCycles=1; % choose how many heart cycles to use
        tIdx = find( Ts.t>tED(cycleStart) & Ts.t<tED(cycleStart+nCycles)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx);
        v=Ts.velocity(tIdx);
        p=Ts.ART(tIdx);%

        %------------- find actual impedance -------------
        V=fft(v);
        P=fft(p);
        Z=P./V;
        N=length(t); dt=t(2)-t(1); fs=1/dt; f=(0:1/N:(1-1/N))*fs;

        %------------- find WindKessel impedance ------------- 
        w=2*pi*f;    % vinkelfrekvens
        R=Z(1);      % should be ~1000
        Adm= 1/R - i*w*C_k; % admitans for paralellkobling av R og C
        Z_WK=1./Adm; % impedans

        % accumulate the error 
        errSum = errSum + abs(abs(Z(2))-abs(Z_WK(2)));
        % have a look at the error of all first 4 frequencies
        err4FreqSum = err4FreqSum + abs(sum(abs(Z(1:4))) - sum(abs(Z_WK(1:4))));
    end
    % add average error to list
    err = [err, errSum/length(cycles)];
    err4Freq = [err4Freq, err4FreqSum/length(cycles)];

end
figure();clf;sgtitle('Compare error of 2. harmonic');
plot(C, err)
xlabel('Capacitance');ylabel('Error')

[min_val, min_idx] = min(err)
actual_compliance = C(min_idx)

figure();clf;sgtitle('Compare error of first 4 frequencies');
plot(C, err4Freq)
xlabel('Capacitance');ylabel('Error')

%%
figure();clf;sgtitle(name);
plot(relPwr);hold on;plot(relPwrWK);hold off;
ylim([0, max([relPwr, relPwrWK]+5)]);