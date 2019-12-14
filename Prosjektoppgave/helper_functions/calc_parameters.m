function [R_arr, C_arr, t_axis] = calc_parameters(ds, C_lim, R_lim)
    %CALC_PARAMETERS Summary of this function goes here
    %   param "limit" is only used if we have poles that explodes in our
    %   impedance and we know what a reasonable limit would be. If above,
    %   set it equal to the previous value.
    if nargin<3 R_lim = log(0); end % infinite
    if nargin<2 C_lim = log(0); end % infinite
    delay = 0.07;
    Ts = ds.Ts;
    tED = ds.t_pulses-delay;
    nBeats = 1; % number of beats we analyse per "cycle"
    totCycles = floor(length(tED/nBeats));
    Tsamp = tED(2)-tED(1); Fsamp = 1/Tsamp;
    C_arr = []; R_arr = []; t_axis = [];
    for k = 1:totCycles-nBeats-1
        tIdx = find( Ts.t>tED(k) & Ts.t<tED(k+nBeats)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx); N=length(t);
        if N<=1 continue; end;
        % Create line between start and end of signal
        v_line = linspace(v(1), v(end), N); v = v - v_line'+mean(v_line);  
        p_line = linspace(p(1), p(end), N); p = p - p_line'+mean(p_line);
        % get parameters
        [C, R] = get_parameters(t, p, v, C_lim, R_lim);
        R_arr = [R_arr, R]; C_arr = [C_arr, C]; t_axis = [t_axis, tED(k)];
    end
end

function [C, R] = get_parameters(t, p, v, C_lim, R_lim)
    Tsamp = t(2)-t(1); fs = 1/Tsamp; N = length(t);
    f = (0:1/(N-1):1)*fs;
    V = fft(v); P = fft(p); Z = P./V;
    R = abs(Z(1)); % zero freq component is the real value of Z
    %if R > R_lim
    C_init = abs( 1i*(1/R-1/Z(2)) / (2*pi*f(2)) );
    % get first 4 frequencies (FF - First Frequencies)
    FF = 1:(1*N/fs+2);
    % Filter first 4 frequencies
    Z_FF = abs(Z(FF)); f_FF = f(FF)';
    % Create function to fit
    Z_eq = @(C_fit, x) 1./abs(1/R + 1i*2*pi*x*C_fit);
    % Fit the function with regards to compliance
    fit_model = fit(f_FF, Z_FF, Z_eq, 'Start',C_init, 'lower',1e-5, 'upper',1e-2);

    C = fit_model.C_fit;
end
