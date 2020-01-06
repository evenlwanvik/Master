function plot_pulses(ds, pause_time)
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
        v_line = linspace(v(1), v(end), N); %v = v - v_line'+mean(v_line);  
        p_line = linspace(p(1), p(end), N); %p = p - p_line'+mean(p_line);
        t_line = linspace(t(1), t(end), N);
        subplot(2,1,1); plot(t,v,t_line,v_line);ylabel('Velocity [cm/s]');xlabel('t [s]')
        subplot(2,1,2); plot(t,p,t_line,p_line);ylabel('Pressure [mmHg]')
        pause(pause_time)
    end
end

