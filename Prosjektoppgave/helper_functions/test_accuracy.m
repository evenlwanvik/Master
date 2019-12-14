%% Script to check the accuracy of the heartrate. 
% Checks if the difference between the start and end point is larger than a
% certain threshold, e.g. if it is more than (max-min)/2 there obiously is
% a large discontinuity

function [velocityMissRate, pressureMissRate] = test_accuracy(ds, divRatio)
    Ts = ds.Ts; Tmean = ds.Tmean; delay=ds.delay; Tmax = ds.Tmax;
    
    totCycles = length(ds.t_pulses);
    numOfVelocityFault = 0;
    numOfPressureFault = 0;
    for k = 1:totCycles-20
        tIdx = find( Ts.t>ds.t_pulses(k) & Ts.t<ds.t_pulses(k+1)); % all samples within tED1 - tED2 window
        t=Ts.t(tIdx); v=Ts.velocity(tIdx); p=Ts.ART(tIdx); N=length(t);
        if N==0 continue; end;
        % Create line between start and end of signal
        velocityThreshold = abs(v(1)-v(end));
        vMax = max(v); vMin = min(v); minMaxDiff = (vMax-vMin)/divRatio;
        if velocityThreshold >= minMaxDiff
            numOfVelocityFault = numOfVelocityFault + 1;
        end
        pressureThreshold = abs(p(1)-p(end));
        pMax = max(p); pMin = min(p); minMaxDiff = (pMax-pMin)/divRatio;
        if pressureThreshold >= minMaxDiff
            numOfPressureFault = numOfPressureFault + 1;
        end
    end    
    velocityMissRate = (numOfVelocityFault/totCycles)*100;
    pressureMissRate = (numOfPressureFault/totCycles)*100;
end