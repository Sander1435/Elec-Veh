% ==============
% TU/e QSS Toolbox
% Vehicle body block initialization script
% ==============

%% Global variables
global h                    % Stepsize [s] from block "Driving Cycle"

%% Constants
g           = 9.81;         % Gravitation constant
rho         = 1.18;         % Air density

%% Parameters conversion
mt2m_f      = mt2m_f/100;	% The user made the input in [%]
r_wheel     = d_wheel/2;	% The user specified the diameter, not the radius
 
%% Driveline selection regen fraction

if drivetype == 1
    f_r = 0.6;
elseif drivetype == 2
    f_r = 0.4;
elseif drivetype == 3
    f_r = 0.9;
end