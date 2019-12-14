function plot_dfts(ds, pause_time)
    %PLOT_PULSES Plot velocity and pressure data every "pause_time" second
    %   ds (dataset) is an object with all necessary data
    figure(2)
    Ts=ds.Ts;
    delay=ds.delay 
    tED=ds.t_pulses-delay;
    N = length(tED)
    for i=1:N
        tIdx = find( Ts.t>tED(i) & Ts.t<tED(i+1)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx); N=length(t);
        v_line = linspace(v(1), v(end), N); v = v - v_line'+mean(v_line);  
        p_line = linspace(p(1), p(end), N); p = p - p_line'+mean(p_line);
        Tsamp = t(2)-t(1); fs = 1/Tsamp; f = (0:1/(N-1):1)*fs;
        V = fft(v); P = fft(p); Z = P./V;
        subplot(2,1,1);
        yyaxis left
        plot(f, abs(P))
        yyaxis right
        plot(f, abs(V));grid();xlim([0, 15])
        subplot(2,1,2);
        plot(f, abs(Z));grid();xlim([0, 15])
        pause(pause_time)
    end
end

