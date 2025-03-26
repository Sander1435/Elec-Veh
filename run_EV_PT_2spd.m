
%% 4AUB10 - Electric and hybrid powertrain design
%  Academic year 2024-2025, Q3
%  Runner for 'EV_2spd.slx' model
clear all; close all; clc;
%warning('off','all');

%% add QSS toolbox to the path
p = genpath('F:\Github\Elec-Veh');
addpath(p);

%% Global variables
global N_sim
global g1
global g2
global scale_EM
global Np
global Ns
global s1
global h
%% Parameters definition
% Model parameters for EV with -spd
x = [11.87 4.42 0.514 30 285 2 2 1 1 57.61];  % set of input design parameters, x
g1              = x(1);       % gear ratio [-]
g2              = x(2);       % gear ratio [-]
scale_EM        = x(3);       % EM scale [-]
Np              = x(4);       % Number of cells connected in parallel [-]
Ns              = x(5);       % Number of cells connected in series [-]
motortype       = x(6);       % Motortype 1 is a PMSM and Motortype 2 is a IM [-] 
motornumberinit = x(7);       % This parameter changes between 2 powers of an engine
BatterySort     = x(8);       % Decides between High Power(2) or Energy cells(1) [-]
CellSort        = x(9);       % Decides the exact cells used for HE from 1 to 2 and for HP 1 to 4
s1              = x(10);      % Shift point of the car
%% Model name
sys_name = 'EV_2spd_4AUB10_2023b';
open(sys_name);

%% Setting up parameters of different componenets

% Gear box 
e_gb    = 0.96;   % internal efficiency [-]
Ploss   = 100;    % stationary losses [W]
wem_min = 1;      % Minimum wheel speed beyond which losses are generated [rad/s]

% Battery 
init_SoC=90;
path= [sys_name,'/Battery (ECM)'];
Battery_type=(["High Energy", "High Power"]);                           
CellTypeHE = (["NMC","NMC_M35A"]);
CellTypeHP = (["LFP","LTO","NCA","NCAVTC"]);
if (BatterySort==1)
    set_param(path,'BattType',Battery_type(BatterySort),'BattChem',CellTypeHE(CellSort));
elseif (BatterySort==2)
    set_param(path,'BattType',Battery_type(BatterySort),'BattChem',CellTypeHP(CellSort));
end
battery_weight=str2num(get_param(path,'bt_weight'));
Battery_size=str2num(get_param(path,'bt_cap'));


% Electric machine
J_em = 0.1;       % rotating inertia [kgm2]
Paux = 0;         % auxilary power losses [W]
path= [sys_name,'/Electric motor 1'];
MotorTypes=(["Permanent Magnet (PMSM)","Induction (AC)"]);
motornumberPM=(["32 kW","58 kW"]);
motornumberAC=(["75 kW","187 kW"]);
if (motortype==1)
    set_param(path,"motortype",MotorTypes(motortype),"motornumberPM",motornumberPM(motornumberinit));
    Pemmax=get_param(path,'MotornumberPM');
    Pemmax=str2num(Pemmax(1:(end-3)))*1e3*scale_EM;
    load("motorstringPM.mat");
elseif (motortype==2)
    set_param(path,"motortype",MotorTypes(motortype),"motornumberAC",motornumberAC(motornumberinit));
    Pemmax=get_param(path,'motornumberAC');
    Pemmax=str2num(Pemmax(1:(end-3)))*1e3*scale_EM;
    load("motorstringAC.mat");
end
load([string(components(motornumberinit))+"at"])
T_EM_max=T_EM_max*scale_EM;
T_EM_col=T_EM_col*scale_EM;
clear components;


% vehicle mass
mem=Pemmax/1e3/1.4;                 % E-machine mass [kg]
mb     = battery_weight;          % battery mass [kg]
mtr    = 45;                                  % 1spd transmission mass
m0     = 1180;                                % curb mass [kg]
mcargo = 0;                                   % cargo mass [kg]
mv     = m0 + mb + mem + mtr + mcargo;        % vehicle mass [kg]

% Vehicle parameters
f_r    = 0.4;       %reg. brake fraction, RWD [-]
lambda = 1.05;      % relative rotating inertia parameter [-]
cr     = 0.0174;     % rolling resistance coefficient [-]
cd     = 0.29;      % air drag coefficient [-]
Af     = 2.38;      % frontal area [m2]
dw     = 0.6996;    % wheel diameter [m]

%% Performance calculations
% Acceleration
g = 9.81;        % Gravitation constant
rho = 1.18;      % Air density
kR = 0.55;       % weight distribution during acceleration
muP = 0.90;      % max. tyre-road adhesion coefficient
Ttmax = kR*muP*mv*g*dw/2;  % max. slip / tractive torque
Ptmax = Pemmax*0.97;
vb = Ptmax./(Ttmax./(dw/2)); % base vehicle speed
vf = 100/3.6;                % final speed acceleration
ta = (lambda*mv*(vb^2+vf^2))./(2*(Ptmax - (2/3)*mv*g*cr*vf - (1/5)*rho*cd*Af*vf^3));

% Max. speed vehicle without overrevving machine
vmax1 = (dw/2)* max(w_EM_max)/g1; % [m/s]
vmax2 = (2*Ptmax/(cd*Af*rho))^(1/3);
if (vmax1 <= vmax2)
    vmax = vmax1;
end
if (vmax2 < vmax1 )
    vmax = vmax2;
end
%% Simulation
results = sim(sys_name);
%close_system(sys_name,0);
% Consider last value of the computed fuel consumption vector
cons_result = results.Consumption(end);

% Check whether cycle could be finished exactly in N_sim computational steps;
% if cycle duration is less than N_sim, set fuel consumption to infinite
if (max(size(results.t)) < N_sim)
    cons_result = Inf;
end

%% Data display
disp(['Energy consumption  = ', num2str(cons_result),' [kWh/100km]']);
disp(['E-machine size      = ', num2str(Pemmax/1e3), ' [kW]']);
disp(['Energy battery pack = ', num2str(Battery_size), ' [kWh]']);
disp(['E-machine mass      = ', num2str(mem), ' [kg]']);
disp(['Battery mass        = ', num2str(mb), ' [kg]']);
disp(['Transmission mass   = ', num2str(mtr), ' [kg]']);
disp(['Vehicle mass        = ', num2str(mv), ' [kg]']);
disp(['Acceleration time   = ', num2str(ta), ' [s]']);
disp(['Max. speed          = ', num2str(vmax*3.6), ' [km/h]']);

%%
switch BatterySort
    case 1  % High Energy
        switch CellSort
            case 1 
                V_cell_max = 4.2; 
            case 2  
                V_cell_max = 4.2;
            otherwise
                V_cell_max = 4.2; 
        end
    case 2  % High Power
        switch CellSort
            case 1  
                V_cell_max = 3.65;
            case 2  
                V_cell_max = 2.7;
            case 3  
                V_cell_max = 4.2;
            case 4  
                V_cell_max = 4.2;
            otherwise
                V_cell_max = 4.2; 
        end
    otherwise
        % If unknown, pick a default or raise a warning
        V_cell_max = 4.2;
end

% 2) Compute the maximum battery pack voltage:
battery_pack_max_voltage = Ns * V_cell_max;

% 3) Display or log it:
disp(['Battery pack max voltage: ', num2str(battery_pack_max_voltage), ' V']);

% 4) Check if it exceeds 1200 V:
if battery_pack_max_voltage > 1200
    warning(['Battery pack max voltage (', ...
             num2str(battery_pack_max_voltage), ' V) exceeds 1200 V limit!']);
end
%% --- Price and TCO Calculations ---
battery_cost_per_kWh = 200;    % [euro/kWh]
motor_cost_per_kW    = 16;      % [euro/kW] for the E-machine
inverter_cost_per_kW = 15;      % [euro/kW] for the inverter
base_vehicle_cost    = 17000;   % [euro] (for a car)

Eb = Battery_size;         % Battery pack size [kWh]
Pm = Pemmax / 1e3;             % E-machine power [kW]
Cb = battery_cost_per_kWh * Eb;      % Battery cost
Cm = motor_cost_per_kW * Pm;           % Motor cost
Ci = inverter_cost_per_kW * Pm;        % Inverter cost
Cv = Cb + Cm + Ci + base_vehicle_cost; % Total component cost

% Operational cost assumptions:
electricity_cost = 0.317;       % [euro/kWh current price]
Dc = 23.25;                    % [km]
Es = cons_result/100 * Dc;     % [kWh/cycle]

% Yearly mileage and economic life
Dy = 20000;   % [km/year]
y  = 5;       % [years]

% Operational (running) cost:
Ce = electricity_cost * Es * (Dy * y / Dc);

% Total Cost-of-Ownership (TCO):
TCO = 0.5 * Cv + Ce;

% Display cost results
disp('--- Cost Calculations ---');
fprintf('Total cost-of-ownership (TCO): %.2f euro\n', TCO);

%% Calculate time step vector
dt = [diff(results.t); 0]; 

% Compute total energy over the simulation cycle:
energy_output = sum(results.P_wheel .* dt);         % Energy delivered at wheels
energy_input  = sum((results.P_EM1 + results.Losses) .* dt);  % Energy drawn from the battery

% Compute overall powertrain efficiency
powertrain_efficiency = energy_output / energy_input;

fprintf('Overall powertrain efficiency: %.2f%%\n', powertrain_efficiency * 100);

%%
% efficiency and average power 
P_wheel = results.P_wheel/1e3;
P_bat = (results.P_EM1+results.Losses)/1e3;

E_wheel = trapz(results.t, P_wheel);
E_bat = trapz(results.t, P_bat);

eff = E_wheel/E_bat;

avg_drive_power = E_wheel/(results.t(end)); % in kW
fprintf('Average drive power: %.2f kW\n', avg_drive_power);
%% --- Range Calculation ---
vehicle_range = Battery_size * 100 / cons_result;  % in km
fprintf('Estimated Vehicle Range: %.2f km\n', vehicle_range);
%% Plots (Simulink model data)
fig = figure;
set(fig,'NumberTitle', 'off')
set(fig,'Name', 'Vehicle speed, torque, power')
subplot(221), plotyy(results.t, results.w_wheel, results.t, 3.6*results.w_wheel*dw/2);
ylabel('\omega_w [rad/s], or v [km/h]'); 
grid;
subplot(222), plot(results.t, results.T_wheel); 
ylabel('T_w [Nm]'); 
grid;
subplot(223), plot(results.t, results.P_wheel/1e3); 
ylabel('P_w [kW]'); 
grid; 
xlabel('Time [s]');
subplot(224), plot(results.t, results.P_wheel/1e3, results.t, results.P_MGB/1e3,...
results.t, results.P_EM1/1e3, results.t, (results.P_EM1+results.Losses)/1e3); grid;
xlabel('Time [s]');
ylabel('[kW]');
legend('P_{wheel}','P_{trans}', 'P_{EM}', 'P_{BT,output}');

fig = figure;
set(fig,'NumberTitle', 'off')
set(fig,'Name', 'Powers')
plot(results.P_wheel/1e3, results.P_MGB/1e3,'.', results.P_wheel/1e3,...
results.P_EM1/1e3,'o', results.P_wheel/1e3, (results.P_EM1+results.Losses)/1e3,'.');
grid; 
legend('P_{trans}', 'P_{EM}', 'P_{BT,output}');
xlabel('P_{wheel} [kW]'); ylabel('Powers [kW]');
%%
fig = figure;
set(fig,'NumberTitle', 'off')
set(fig,'Name', 'Battery state of charge')
plot(linspace(min(results.t),max(results.t),length(results.SoC)), results.SoC);
grid;
ylabel('SoC [%]');
xlabel('Time [s]');
ylim([0 100]);

%%
set(fig,'NumberTitle', 'off')
set(fig,'Name', 'Operation points')
v = 0:0.05:2;
[x,y]=meshgrid(w_EM_row*30/pi, T_EM_col);
[c,h] = contourf(x,y, eta_EM_mapM',v);
clabel(c,h);
hold on;
plot(w_EM_max*30/pi, T_EM_max, w_EM_max*30/pi, -T_EM_max, 'r', 'linewidth', 2);
plot(results.w_EM1*30/pi, results.T_EM1, 'o');
xlabel('w_{EM} [rpm]'); ylabel('T_{EM} [Nm]');
grid;

%%
% Motor Speed vs. Time

gear1_idx = results.w_EM1 < s1;  % operating in gear 1
gear2_idx = results.w_EM1 >= s1; % operating in gear 2

figure;
hold on;
plot(results.t(gear1_idx), results.w_EM1(gear1_idx), 'bo', 'DisplayName', 'Gear 1');
plot(results.t(gear2_idx), results.w_EM1(gear2_idx), 'rs', 'DisplayName', 'Gear 2');
xlabel('Time [s]');
ylabel('Motor Speed [rad/s]');
title('Motor Speed vs. Time (Gear 1 vs. Gear 2)');
legend('Location','best');
grid on;
hold off;

%% Time Analysis per Gear ---

dt = [diff(results.t); 0];

gear1_mask = results.w_EM1 < s1;
gear2_mask = results.w_EM1 >= s1;

% Total time spent in each gear
timeGear1 = sum(dt(gear1_mask));
timeGear2 = sum(dt(gear2_mask));

fprintf('Time in Gear 1: %.2f s\n', timeGear1);
fprintf('Time in Gear 2: %.2f s\n', timeGear2);

