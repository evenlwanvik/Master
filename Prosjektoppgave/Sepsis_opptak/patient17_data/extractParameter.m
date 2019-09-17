%% Create function to extract parameters

%load 20190117T145728_IQ_Sepsis-4min_traces;
load 20190123T105641_IQ_Sepsis-4min_traces;

%figure();clf;
%hold on
[R_arr, C_arr, t_axis] = loop_through_dataset(Ts, Tmean, 3);
%hold off
% Average?
%R_arr = movmean(R_arr, 3); C_arr = movmean(C_arr, 3);
figure();clf;
yyaxis left
plot(t_axis, R_arr)
yyaxis right
plot(t_axis, C_arr)
    
function [R_arr, C_arr, t_axis] = loop_through_dataset(Ts, Tmean, nCycles)
    delay = 0.07;
    tED = Tmean.tED-delay;
    totCycles = length(tED)
    
    C_arr = []; R_arr = []; t_axis = [];

    k = 1;
    while (k+nCycles) < totCycles
        tIdx = find( Ts.t>tED(k) & Ts.t<tED(k+nCycles)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); 
        v=Ts.velocity(tIdx);
        p=Ts.ART(tIdx);
        [C, R] = get_parameters(t, p, v);
        R_arr = [R_arr, R];
        C_arr = [C_arr, C];
        t_axis = [t_axis, tED(k)];
        k = k + nCycles;
    end
end

function [C, R] = get_parameters(t, p, v)

    Tsamp = t(2)-t(1); fs = 1/Tsamp;
    N = length(t);
    f = (0:1/(N-1):1)*fs;
    V = fft(v); 
    P = fft(p); 
    Z = P./V;
    R = Z(1); % zero freq component is the real value of Z

    % find the first harmonic frequency of our signal
    [x, i_Fv] = max(abs(V(2:N-1))); Fv = f(i_Fv-1);
    [x, i_Fp] = max(abs(P(2:N-1))); Fp = f(i_Fp-1);
    %[x, i_Fz] = max(abs(Z(2:N-1))); Fz = f(i_Fz-1);

    %C = abs(Z(i_Wz)/(1i*Fz));
    C = abs( ( 1/R - 1/Z(i_Fp) ) / (Fp*2*pi) )
    
    %Adm = 1/R + i*Wz*C;
    %Q=Adm.*P;
    %q=abs(real(ifft(Q)));% flow i tidsplan
    %subplot(1,3,1);plot(t,p);title('Heart cycle - pressure')
    %subplot(1,3,2);plot(t,v);title('Heart cycle - velocity')
    %subplot(1,3,3);plot(t,q);title('Velocity remade')
end

