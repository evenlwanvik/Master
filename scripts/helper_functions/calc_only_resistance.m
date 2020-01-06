function [R_arr, t_axis] = calc_only_resistance(ds)
    %CALC_PARAMETERS Summary of this function goes here
    %   Same as calc_parameters but we only calculate resistance for more than one pulse.
    if nargin<3 R_lim = log(0); end % infinite
    if nargin<2 C_lim = log(0); end % infinite
    delay = 0.07;
    Ts = ds.Ts;
    tED = ds.t_pulses-delay;
    nBeats = 3; % number of beats we analyse per "cycle"
    totCycles = floor(length(tED/nBeats));
    R_arr = []; t_axis = [];
    for k = 1:nBeats:totCycles-nBeats-1
        tIdx = find( Ts.t>tED(k) & Ts.t<tED(k+nBeats)); % all samples within tED1 - tED2 window
        pause(0.2)
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx); N=length(t);
        if N<=1 continue; end;
        % Create line between start and end of signal
        v_line = linspace(v(1), v(end), N); v = v - v_line'+mean(v_line);  
        p_line = linspace(p(1), p(end), N); p = p - p_line'+mean(p_line);
        % get parameters
        R = get_parameters(t, p, v);
        R_arr = [R_arr, R]; t_axis = [t_axis, tED(k)];
    end
end

function [R] = get_parameters(t, p, v)
    Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);
    V = fft(v); P = fft(p); Z = P./V;
    R = abs(Z(1)); % zero freq component is the real value of Z
end