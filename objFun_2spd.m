function [fVal, overspeed] = objFun_2spd(x)
% OBJFUN_2SPD  Objective function for optimizing run_EV_PT_2spd.m
%
% This version puts *all* variables used by the blocks into the base
% workspace before calling the base script. Then the blocks that reference
% them (like init_SoC, Np, Ns, scale_EM, etc.) can find them.

    % By default we won't detect overspeed unless an error is thrown
    overspeed = false;
    fVal = 1e9;

    % If you want to vary 4 parameters: [g1, g2, scaleEM, shiftPoint],
    % let's say x has exactly 4 elements:
    if length(x) < 4
        error('objFun_2spd: input x must have >= 4 elements: [g1, g2, scale_EM, shift].');
    end

    %% 1) Extract your design variables from x
    %   We'll let these four be the ones we optimize.
    g1_val     = x(1);
    g2_val     = x(2);
    scale_val  = x(3);
    s1_val     = x(4);

    %% 2) Provide default / fixed values for the other parameters
    %   The model references them by name (init_SoC, Np, Ns, etc.)
    %   If you want to optimize them too, just read them from x() also.
    init_SoC = 90;            % e.g. SoC start
    Np_val   = 16;            % #parallel
    Ns_val   = 96;            % #series
    BatterySort_val = 1;      % 1=HighEnergy
    CellSort_val    = 1;      % e.g. NMC
    motortype_val   = 1;      % 1=PMSM
    motornumberinit_val = 1;  % which motor rating
    Ploss_val       = 100;    % for the gearbox
    wem_min_val     = 1;
    e_gb_val        = 0.97;

    % If needed, you can also define any other variables the blocks reference
    % e.g. "mem", "mv", etc. Usually, you compute them in the base script,
    % so we only need to pass in the "inputs" to that script. The script
    % can then compute mem, mv, etc. But if the blocks themselves want e.g. 'mv'
    % directly, you'd also do assignin('base','mv', someValue).

    %% 4) Run the base script inside a specialized try/catch
    %    (only handle overspeed; rethrow everything else)
    try
        run('run_EV_PT_2spd.m');
    catch ME
        if strcmp(ME.identifier,'Simulink:Engine:OverSpeed')
            fVal = Inf;
            overspeed = true;
            return;
        else
            rethrow(ME); % rethrow other errors to see them
        end
    end

    %% 5) Retrieve final consumption from base workspace
    if evalin('base','exist(''cons_result'',''var'')')
        fVal = evalin('base','cons_result');
    else
        fVal = Inf;
        return;
    end

    %% 6) Check constraints
    %    We assume top speed is in m/s, so multiply by 3.6 for the 150 km/h constraint.
    ta_req   = 11;   % must do 0-100 in <= 11s
    vmax_req = 150;  % must be >= 150 km/h

    if evalin('base','exist(''ta'',''var'')')
        ta_sim = evalin('base','ta');
    else
        ta_sim = 999;
    end

    if evalin('base','exist(''vmax'',''var'')')
        vmax_mps = evalin('base','vmax');
        vmax_kmh = vmax_mps * 3.6;
    else
        vmax_kmh = 0;
    end

    if ta_sim > ta_req
        fVal = Inf;
    end

    if vmax_kmh < vmax_req
        fVal = Inf;
    end

end
