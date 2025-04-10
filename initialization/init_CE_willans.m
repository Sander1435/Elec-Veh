% TU/e QSS Toolbox
% Combustion Engine block (based on Willans-approximation) initialization script
% ==============

%% Global variables
global h                % Stepsize [s] from block "Driving Cycle"

%% Torque at idle
T_CE_idle = P_CE_idle / w_CE_idle;      % [Nm]

%% Engine
V_d     = V_d/1000;                     % [l] -> [m^3]
S       = S/1000;                       % [mm] -> [m]
B       = B/1000;                       % [mm] -> [m]

%% Which engine type?
switch engine_type
    case 1                                  % -> Otto      
        k1  = 1.44 * 1e5;                   % [Pa]
        k2  = 0.46;                         % [-]
        k3  = 9.1 * 1e-4;                   % [s^2/m^2]
        k4  = 0.075;                        % [m]
    case 2                                  % -> Diesel      
        k1  = 1.44 * 1e5;                   % [Pa]
        k2  = 0.50;                         % [-]
        k3  = 1.1 * 1e-3;                   % [s^2/m^2]
        k4  = 0.075;                        % [m]
end
