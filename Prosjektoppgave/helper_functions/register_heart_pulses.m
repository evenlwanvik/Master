function [t_pulses] = register_heart_pulses(ecg, time, limit)
    %REGISTER_HEART_PULSES If heart pulses are wrongfully registered, use
    %this function to create our own Tmean time axis
    %   If ecg goes above a certain "limit", we have a pulse. You should
    %   investigate what the limit should be by plotting the full ceg arr.
    t_pulses = [];
    above=false;
    for i=1:length(time)
        if ecg(i) > limit
            if above==false
                t_pulses = [t_pulses, time(i)];
                above=true;
            end
        else
            above=false;
        end
    end
end

