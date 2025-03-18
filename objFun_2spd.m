function [fVal, overspeed] = objFun_2spd(x)
% OBJFUN_2SPD  Objective function for optimizing run_EV_PT_2spd.m
%
% x = [g1, g2, scale_EM, shift]
% Objective: minimize consumption [kWh/100km]
% Constraints: top speed >= 150 km/h, acceleration time <= 11 s.
% If the simulation throws an OverSpeed or Overload error, or does not complete, return Inf.

    overspeed = false;
    fVal = 1e9;

    if length(x) < 4
        error('objFun_2spd: input x must have >= 4 elements: [g1, g2, scale_EM, shift].');
    end

    % Clear previous base workspace variables (except global N_sim)
    evalin('base','clearvars -except N_sim');

    % Extract design variables
    g1_val = x(1);
    g2_val = x(2);
    scale_val = x(3);
    s1_val = x(4);

    % Build base vector and assign to base workspace
    baseX = [g1_val, g2_val, scale_val, s1_val];
    assignin('base','x', baseX);

    % Run simulation and catch OverSpeed/Overload errors
    try
        run('run_EV_PT_2spd.m');
    catch ME
        if contains(ME.message, 'OverSpeed') || contains(ME.message, 'Overload')
            fVal = Inf;
            overspeed = true;
            return;
        else
            rethrow(ME);
        end
    end

    % Retrieve consumption (kWh/100km)
    if evalin('base','exist(''cons_result'',''var'')')
        consVal = evalin('base','cons_result');
    else
        fVal = Inf;
        return;
    end

    % Check that the simulation completed the entire cycle
    if evalin('base','exist(''results'',''var'')')
        simTime = evalin('base','results.t');
        global N_sim
        if length(simTime) < N_sim
            fVal = Inf;
            return;
        end
    end

    % Retrieve acceleration time and top speed
    if evalin('base','exist(''ta'',''var'')')
        ta_val = evalin('base','ta');
    else
        ta_val = 999;
    end
    if evalin('base','exist(''vmax'',''var'')')
        vmax_mps = evalin('base','vmax');
        vmax_kmh = vmax_mps * 3.6;
    else
        vmax_kmh = 0;
    end

    % Constraint requirements
    ta_req = 11;    % seconds
    vmax_req = 150; % km/h

    if ta_val > ta_req || vmax_kmh < vmax_req
        fVal = Inf;
        return;
    end

    % If all constraints are met, objective is consumption.
    fVal = consVal;
end
