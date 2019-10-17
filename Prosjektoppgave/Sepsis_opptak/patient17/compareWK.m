%% Compare the first 5 Hz of actual impedance and WK model impedance
addpath('Master/Prosjektoppgave/Sepsis_opptak/patient17/data') 
%load 20190117T145728_IQ_Sepsis-4min_traces;
%date = '17.01.2019'
load 20190123T105641_IQ_Sepsis-4min_traces;
date = '23.01.2019'
name = append('Sepsis 4min_traces ', date)

relPwr = []; relPwrWK = []; 

for cycleStart=1:100:300
    figure();clf;sgtitle(name);
    %-------- measured data from single heart cycle --------

    delay = 0.07;
    tED = Tmean.tED-delay;
    %cycleStart=33; % choose heart cycle (~300 possibilities)
    nCycles=1; % choose how many heart cycles to use
    tIdx = find( Ts.t>tED(cycleStart) & Ts.t<tED(cycleStart+nCycles)); % all samples within tED1 - tED2 window
    t=Ts.t(tIdx);
    v=Ts.velocity(tIdx);
    p=Ts.ART(tIdx);

    %------------- find actual impedance -------------

    V=fft(v);
    P=fft(p);
    Z=P./V;
    N=length(t); dt=t(2)-t(1); fs=1/dt; f=(0:1/N:(1-1/N))*fs;

    %------------- find WindKessel impedance ------------- 
    
    w=2*pi*f;    % vinkelfrekvens
    R=Z(1);      % should be ~1000
    C = abs( i*(1/Z(2)-1/R) / (2*pi*f(2)) ); % should be ~1.18e-4
    C=0.0002;
    Adm= 1/R - i*w*C; % admitans for paralellkobling av R og C
    Z_WK=1./Adm; % impedans

    plot(f,abs(Z),f,abs(Z_WK));
    %xlim([0,10]) 

    pwrSum = sum(abs(Z(1:4)));
    relPwr = [abs(Z(2))*100/pwrSum, relPwr]
    pwrSumWK = sum(abs(Z_WK(1:4)));
    relPwrWK = [abs(Z_WK(2))*100/pwrSumWK, relPwrWK]
end
    
figure();clf;sgtitle(name);
plot(relPwr);hold on;plot(relPwrWK);hold off;
%ylim([0, ]);%max([relPwr, relPwrWK]+5)]);
