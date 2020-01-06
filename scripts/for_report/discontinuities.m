figure(18);     
tED = Tmean3.tED-delay3;
t = Ts3.t; p = Ts3.ART; v = Ts3.velocity;
tIdx = find( t>tED(200) & t<tED(200+1));

p = p(tIdx); v = v(tIdx); t = t(tIdx);
N = length(t); fs = 1/(t(2)-t(1)); f = (0:1/(N-1):1)*fs;

v_line = linspace(v(1), v(end), N); v_new = v - v_line'+mean(v_line);  
p_line = linspace(p(1), p(end), N); p_new = p - p_line'+mean(p_line);


%%
t_line = linspace(t(1), t(end), N);
straight = linspace(p(1), p(1), N);
plot(t,p,t_line,p_line, '--',t_line,straight, 'k');ylabel('Pressure [mmHg]');xlabel('Time [sec]')

%%

P = fft(p); V = fft(v); Z = P./V/N;
P = fft(p_new); V = fft(v_new); Z_new = P./V/N;
plot(f, abs(Z), f, (Z_new))
xlim([0,20])
xlabel('Frequency [Hz]')
ylabel('|Z(\omega)|')
legend('w/ compensated')
legend('FFT of unchanged signal', 'FFT of continious signal')