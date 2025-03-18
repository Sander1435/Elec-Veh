%% 4AUB10 - Electric and Hybrid Powertrain Design
% Academic year 2024-2025, Q3
% Runner for 'EV_2spd.slx' model (2-speed EV powertrain)
close all; clc;
%warning('off','all');

%% Add QSS toolbox to the path
p = genpath('QSS_TB');
addpath(p);

%% Global variables
global N_sim
N_sim = 1800;  % Number of simulation steps

%% Parameters Definition
% Design vector: [g1, g2, scale_EM, shift]
if ~exist('x','var')
    x = [6.5, 3.25, 1.98, 28];  % Default design vector if none provided
end

% Extract design variables
g1       = x(1);   % Gear ratio #1
g2       = x(2);   % Gear ratio #2
scale_EM = x(3);   % Motor scale factor
s1       = x(4);   % Shift threshold

% Fixed battery parameters (optimization limited to drivetrain)
Np = 88;  % Number of cells in parallel
Ns = 96;  % Number of cells in series

% Other fixed parameters
motortype       = 2;   % 1: PMSM, 2: Induction
motornumberinit = 1;   % Base motor index
BatterySort     = 1;   % 1 = High Energy; 2 = High Power
CellSort        = 2;   % Specific cell selection index
init_SoC        = 90;
e_gb            = 0.97;    % Gearbox internal efficiency
Ploss           = 100;     % Stationary losses [W]
wem_min         = 1;       % Minimum wheel speed for loss generation [rad/s]

%% Battery Specifications
cellVoltage  = 3.7;    % [V]
cellCapacity = 3.4;    % [Ah]
cellMass     = 0.045;  % [kg per cell]
battery_size_kWh = (Np * Ns * cellVoltage * cellCapacity) / 1000; 
battery_weight   = (Np * Ns * cellMass);

%% Place Variables into the Base Workspace
assignin('base','x', x);
assignin('base','g1', g1);
assignin('base','g2', g2);
assignin('base','scale_EM', scale_EM);
assignin('base','s1', s1);
assignin('base','Np', Np);
assignin('base','Ns', Ns);
assignin('base','motortype', motortype);
assignin('base','motornumberinit', motornumberinit);
assignin('base','BatterySort', BatterySort);
assignin('base','CellSort', CellSort);
assignin('base','init_SoC', init_SoC);
assignin('base','e_gb', e_gb);
assignin('base','Ploss', Ploss);
assignin('base','wem_min', wem_min);
assignin('base','battery_size_kWh', battery_size_kWh);

%% Model Name and Open Simulink Model
sys_name = 'EV_2spd_4AUB10_2023b';
open_system(sys_name);

%% Setting Up Battery Block Parameters
pathBatt = [sys_name, '/Battery (ECM)'];
Battery_type = (["High Energy", "High Power"]);
CellTypeHE   = (["NMC", "NMC_M35A"]);
CellTypeHP   = (["LFP", "LTO", "NCA", "NCAVTC"]);
if (BatterySort == 1)
    set_param(pathBatt, 'BattType', Battery_type(BatterySort), 'BattChem', CellTypeHE(CellSort));
else
    set_param(pathBatt, 'BattType', Battery_type(BatterySort), 'BattChem', CellTypeHP(CellSort));
end

%% Motor Block Setup
pathMotor = [sys_name, '/Electric motor 1'];
MotorTypes = (["Permanent Magnet (PMSM)", "Induction (AC)"]);
motornumberPM = (["32 kW", "58 kW"]);
motornumberAC = (["75 kW", "187 kW"]);
if (motortype == 1)
    set_param(pathMotor, "motortype", MotorTypes(motortype), "motornumberPM", motornumberPM(motornumberinit));
    Pemmax = get_param(pathMotor, 'MotornumberPM');
    Pemmax = str2num(Pemmax(1:(end-3)))*1e3*scale_EM;
    load("motorstringPM.mat");
else
    set_param(pathMotor, "motortype", MotorTypes(motortype), "motornumberAC", motornumberAC(motornumberinit));
    Pemmax = get_param(pathMotor, 'motornumberAC');
    Pemmax = str2num(Pemmax(1:(end-3)))*1e3*scale_EM;
    load("motorstringAC.mat");
end
load([string(components(motornumberinit))+"at"])
T_EM_max = T_EM_max * scale_EM;
T_EM_col = T_EM_col * scale_EM;
clear components;

%% Vehicle Mass Calculations
mem   = Pemmax / 1e3 / 1.4;  % E-machine mass [kg]
mtr   = 17;                % Transmission mass [kg]
m0    = 1180;              % Curb mass [kg]
mcargo = 0;                % Cargo mass [kg]
mb    = battery_weight;    % Battery mass [kg]
mv    = m0 + mb + mem + mtr + mcargo;  % Total vehicle mass [kg]
assignin('base','mem', mem);
assignin('base','mb', mb);
assignin('base','mtr', mtr);
assignin('base','m0', m0);
assignin('base','mcargo', mcargo);
assignin('base','mv', mv);

%% Rolling and Aerodynamic Parameters
f_r    = 0.4;
lambda = 1.05;
cr     = 0.0174;
cd     = 0.29;
Af     = 2.38;
dw     = 0.6996;
assignin('base','f_r', f_r);
assignin('base','lambda', lambda);
assignin('base','cr', cr);
assignin('base','cd', cd);
assignin('base','Af', Af);
assignin('base','dw', dw);

%% Performance Calculations
g   = 9.81;        % Gravitational constant [m/s^2]
rho = 1.18;        % Air density [kg/m^3]
kR  = 0.55;        % Weight distribution factor during acceleration
muP = 0.90;        % Maximum tyre-road adhesion coefficient
Ttmax = kR * muP * mv * g * dw / 2;  % Maximum tractive torque [Nm]
Ptmax = Pemmax * 0.97;              % Maximum power [W]
vb    = Ptmax / (Ttmax / (dw/2));    % Base vehicle speed [m/s]
vf    = 100/3.6;                     % Final speed for acceleration [m/s]
ta    = (lambda * mv * (vb^2 + vf^2)) / (2 * (Ptmax - (2/3) * mv * g * cr * vf - (1/5) * rho * cd * Af * vf^3));
assignin('base','g', g);
assignin('base','rho', rho);
assignin('base','kR', kR);
assignin('base','muP', muP);
assignin('base','Ttmax', Ttmax);
assignin('base','Ptmax', Ptmax);
assignin('base','vb', vb);
assignin('base','vf', vf);
assignin('base','ta', ta);

% Maximum speed using 2nd gear
vmax = (dw/2) * max(w_EM_max) / g2;
assignin('base','vmax', vmax);

%% Simulation
results = sim(sys_name);
cons_result = results.Consumption(end);  % Energy consumption [kWh/100km]

%% Data Display
disp(['Energy consumption  = ', num2str(cons_result), ' [kWh/100km]']);
disp(['E-machine size      = ', num2str(Pemmax/1e3), ' [kW]']);
disp(['Battery size        = ', num2str(battery_size_kWh), ' [kWh]']);
disp(['E-machine mass      = ', num2str(mem), ' [kg]']);
disp(['Battery mass        = ', num2str(mb), ' [kg]']);
disp(['Transmission mass   = ', num2str(mtr), ' [kg]']);
disp(['Vehicle mass        = ', num2str(mv), ' [kg]']);
disp(['Acceleration time   = ', num2str(ta), ' [s]']);
disp(['Max. speed          = ', num2str(vmax*3.6), ' [km/h]']);

%% --- Price and TCO Calculations ---
% Assumptions based on the assignment for a compact vehicle (BMW i3)
battery_cost_per_kWh = 200;    % [euro/kWh]
motor_cost_per_kW    = 16;      % [euro/kW] for the E-machine
inverter_cost_per_kW = 15;      % [euro/kW] for the inverter
base_vehicle_cost    = 17000;   % [euro] (for a car)

% Component costs (using battery pack size and motor power)
Eb = battery_size_kWh;         % Battery pack size [kWh]
Pm = Pemmax / 1e3;             % E-machine power [kW]
Cb = battery_cost_per_kWh * Eb;      % Battery cost
Cm = motor_cost_per_kW * Pm;           % Motor cost
Ci = inverter_cost_per_kW * Pm;        % Inverter cost
Cv = Cb + Cm + Ci + base_vehicle_cost; % Total component cost

% Operational cost assumptions:
electricity_cost = 0.50;       % [euro/kWh]
% Assume a drive cycle distance (Dc) as in the WLTP example
Dc = 23.25;                    % [km]
% Energy consumption is given in [kWh/100km] so compute energy used per cycle:
Es = cons_result/100 * Dc;     % [kWh/cycle]

% Yearly mileage and economic life for a BMW i3
Dy = 20000;   % [km/year]
y  = 5;       % [years]

% Operational (running) cost:
Ce = electricity_cost * Es * (Dy * y / Dc);

% Total Cost-of-Ownership (TCO):
% Here, TCO is defined as 50% of the component cost plus the operational cost.
TCO = 0.5 * Cv + Ce;

% Sales price assumptions:
profit_margin = 0.10;            % 10% profit margin
sales_price   = (1 + profit_margin) * TCO;
VAT           = 0.21;            % 21% tax (VAT)
final_price   = sales_price * (1 + VAT);

assignin('base','final_price', final_price);


% Display cost results
disp('--- Cost Calculations ---');
fprintf('Battery cost: %.2f euro\n', Cb);
fprintf('Motor cost: %.2f euro\n', Cm);
fprintf('Inverter cost: %.2f euro\n', Ci);
fprintf('Base vehicle cost: %.2f euro\n', base_vehicle_cost);
fprintf('Total component cost, Cv: %.2f euro\n', Cv);
fprintf('Operational cost, Ce: %.2f euro\n', Ce);
fprintf('Total cost-of-ownership (TCO): %.2f euro\n', TCO);
fprintf('Sales price (with 10%% profit): %.2f euro\n', sales_price);
fprintf('Final price (including 21%% VAT): %.2f euro\n', final_price);

%% --- Range Calculation ---
% Compute the vehicle range based on the battery pack energy and energy consumption.
% cons_result is [kWh/100km] and battery_size_kWh is in kWh.
vehicle_range = battery_size_kWh * 100 / cons_result;  % in km
fprintf('Estimated Vehicle Range: %.2f km\n', vehicle_range);
assignin('base','vehicle_range', vehicle_range);
%% --- Plots ---
fig = figure;
set(fig, 'NumberTitle', 'off', 'Name', 'Vehicle speed, torque, power');
subplot(221), plotyy(results.t, results.w_wheel, results.t, 3.6*results.w_wheel*dw/2);
ylabel('\omega_w [rad/s] or v [km/h]'); grid;
subplot(222), plot(results.t, results.T_wheel); ylabel('T_w [Nm]'); grid;
subplot(223), plot(results.t, results.P_wheel/1e3); ylabel('P_w [kW]'); grid; xlabel('Time [s]');
subplot(224), plot(results.t, results.P_wheel/1e3, results.t, results.P_MGB/1e3, ...
    results.t, results.P_EM1/1e3, results.t, (results.P_EM1+results.Losses)/1e3);
xlabel('Time [s]'); ylabel('[kW]'); grid; 
legend('P_{wheel}','P_{trans}','P_{EM}','P_{BT,output}');

fig = figure;
set(fig, 'NumberTitle', 'off', 'Name', 'Powers');
plot(results.P_wheel/1e3, results.P_MGB/1e3, '.', ...
    results.P_wheel/1e3, results.P_EM1/1e3, 'o', ...
    results.P_wheel/1e3, (results.P_EM1+results.Losses)/1e3, '.');
xlabel('P_{wheel} [kW]'); ylabel('Powers [kW]'); grid; 
legend('P_{trans}', 'P_{EM}', 'P_{BT,output}');

fig = figure;
set(fig, 'NumberTitle', 'off', 'Name', 'Battery state of charge');
plot(linspace(min(results.t), max(results.t), length(results.SoC)), results.SoC);
xlabel('Time [s]'); ylabel('SoC [%]'); grid; ylim([0 100]);

fig = figure;
set(fig, 'NumberTitle', 'off', 'Name', 'Operation points');
v = 0:0.05:2;
[xv, yv] = meshgrid(w_EM_row*30/pi, T_EM_col);
[c, h] = contourf(xv, yv, eta_EM_mapM', v); clabel(c, h); hold on;
plot(w_EM_max*30/pi, T_EM_max, w_EM_max*30/pi, -T_EM_max, 'r', 'linewidth', 2);
plot(results.w_EM1*30/pi, results.T_EM1, 'o');
xlabel('w_{EM} [rpm]'); ylabel('T_{EM} [Nm]'); grid;
