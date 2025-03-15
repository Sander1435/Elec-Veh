function [fVal] = objFun_2spd(x)
% objFun_2spd Objective function for optimizing EV_2spd_4AUB10_2023b
%
% INPUT:
%   x : design vector, e.g. [g1, g2, scale_EM, Np, Ns, motortype, motornumberinit, BatterySort, CellSort, shiftPoint]
%
% OUTPUT:
%   fVal : the "cost" to be minimized (kWh/100km or Inf if constraints not met)

%% Make sure we reference the global variable for N_sim if you use it
global N_sim

% --- 1) EXTRACT DESIGN VARIABLES FROM x
g1              = x(1);   % 1st gear ratio
g2              = x(2);   % 2nd gear ratio
scale_EM        = x(3);   % motor scale factor
Np              = x(4);   % # parallel battery cells
Ns              = x(5);   % # series battery cells
motortype       = x(6);   % 1=PMSM, 2=AC Induction
motornumberinit = x(7);   % index to pick e.g. "32 kW", "58 kW", etc.
BatterySort     = x(8);   % 1=High Energy, 2=High Power
CellSort        = x(9);   % e.g. NMC=1, LFP=1, ...
s1              = x(10);  % shift threshold (motor speed)

% --- 2) SET ANY HARD-CODED REQUIREMENTS
ta_req    = 9;    % e.g. must accelerate 0-100 km/h in under 9s
vmax_req  = 140;  % must reach at least 140 km/h
% You can adjust these to your actual constraints

% --- 3) LOAD / OPEN SIMULINK MODEL
sys_name = 'EV_2spd_4AUB10_2023b';
% Instead of open(...), do "load_system(...)" in an optimization to avoid popups
load_system(sys_name);

% --- 4) PUSH PARAMETERS INTO WORKSPACE (or set them in the script directly)
% We do this by assigning them in the base workspace. Another approach is
% to define them as global variables if your code references them that way.
assignin('base','g1',g1);
assignin('base','g2',g2);
assignin('base','scale_EM',scale_EM);
assignin('base','Np',Np);
assignin('base','Ns',Ns);
assignin('base','motortype',motortype);
assignin('base','motornumberinit',motornumberinit);
assignin('base','BatterySort',BatterySort);
assignin('base','CellSort',CellSort);
assignin('base','s1',s1);

% If your code references other variables (like init_SoC=90), you can do:
%   assignin('base','init_SoC',90);
% Or if your code sets them directly, that’s fine too.

% --- 5) RUN THE SIMULATION
try
    simOut = sim(sys_name, 'SrcWorkspace','base');
catch ME
    % If an error occurs (like overspeed) => cost = Inf
    fVal = Inf;
    return
end

% If the sim ran, retrieve final consumption:
cons_result = simOut.Consumption(end);

% Check if simulation ended before finishing
if length(simOut.t) < N_sim
    fVal = Inf;
    return
end

% --- 6) RETRIEVE PERFORMANCE FIGURES
ta    = simOut.ta(1);    % your script or workspace might store acceleration time in simOut
vmax  = simOut.vmax(1);  % top speed in m/s or km/h, depends on your code

% If your model calculates them differently, read or compute them accordingly,
% e.g. from your "ta = (lambda*mv*(vb^2+vf^2))/(...)" expression in the script.

% Convert vmax to km/h if it’s in m/s
%   vmax_kmh = vmax * 3.6;

% --- 7) CHECK CONSTRAINTS
if ta > ta_req
    fVal = Inf;
    return
end

% If your code’s "vmax" is in m/s, compare to (vmax_req / 3.6). Otherwise if it’s in km/h, compare directly
if vmax < vmax_req
    fVal = Inf;
    return
end

% --- 8) IF FEASIBLE => fVal = consumption
fVal = cons_result;   % e.g. kWh/100 km
