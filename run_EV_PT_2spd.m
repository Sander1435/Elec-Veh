%% 4AUB10 - Electric and hybrid powertrain design
%  Academic year 2024-2025, Q3
%  Runner for 'EV_2spd.slx' model
clear all; close all; clc;
%warning('off','all');

%% add QSS toolbox to the path
p = genpath('C:\Users\smaek\Desktop\TuE\Elec and Hybrid\QSS_TB');
addpath(p);

%% Global variables
global N_sim

%% Parameters definition
% --------------------------------------------------------------
% Adjusted for a higher second gear & earlier shift:
%   g1=12, g2=6, scale_EM=4, Np=16, Ns=96,
%   motortype=1 => PMSM, motornumberinit=1 => "32 kW",
%   BatterySort=1 => High Energy, CellSort=1 => NMC,
%   s1=300 => shift around 300 rad/s motor speed
% --------------------------------------------------------------
x = [12 6 4 16 96 1 1 1 1 300];
g1              = x(1);       % 1st gear ratio
g2              = x(2);       % 2nd gear ratio
scale_EM        = x(3);       % EM power scale factor
Np              = x(4);       % # parallel cells
Ns              = x(5);       % # series cells
motortype       = x(6);       % 1=PMSM, 2=AC induction
motornumberinit = x(7);       % pick "32 kW" base
BatterySort     = x(8);       % 1=High Energy, 2=High Power
CellSort        = x(9);       % e.g. NMC
s1              = x(10);      % shift point (motor speed threshold)

%% Model name
sys_name = 'EV_2spd_4AUB10_2023b';
open(sys_name);

%% Setting up parameters of different components

% Gear box
e_gb    = 0.97;
Ploss   = 100;
wem_min = 1;

% Battery
init_SoC=90;
path= [sys_name,'/Battery (ECM)'];
Battery_type=(["High Energy", "High Power"]);
CellTypeHE = (["NMC","NMC_M35A"]);
CellTypeHP = (["LFP","LTO","NCA","NCAVTC"]);
if (BatterySort==1)
    set_param(path,'BattType',Battery_type(BatterySort),'BattChem',CellTypeHE(CellSort));
else
    set_param(path,'BattType',Battery_type(BatterySort),'BattChem',CellTypeHP(CellSort));
end
battery_weight = str2num(get_param(path,'bt_weight'));
Battery_size   = str2num(get_param(path,'bt_cap'));

% Electric machine
J_em = 0.1;  % inertia
Paux= 0;     % aux losses
path= [sys_name,'/Electric motor 1'];
MotorTypes=(["Permanent Magnet (PMSM)","Induction (AC)"]);
motornumberPM=(["32 kW","58 kW"]);
motornumberAC=(["75 kW","187 kW"]);

if (motortype==1)
    set_param(path,"motortype",MotorTypes(motortype),"motornumberPM",motornumberPM(motornumberinit));
    Pemmax = get_param(path,'motornumberPM');
    Pemmax = str2num(Pemmax(1:(end-3)))*1e3*scale_EM;   % e.g. 32*4=128 kW
    load("motorstringPM.mat");
else
    set_param(path,"motortype",MotorTypes(motortype),"motornumberAC",motornumberAC(motornumberinit));
    Pemmax = get_param(path,'motornumberAC');
    Pemmax = str2num(Pemmax(1:(end-3)))*1e3*scale_EM;
    load("motorstringAC.mat");
end

load([string(components(motornumberinit))+"at"])  % loads T_EM_max, T_EM_col
T_EM_max = T_EM_max*scale_EM;
T_EM_col = T_EM_col*scale_EM;
clear components;

% Vehicle mass (BMW i3)
mem   = Pemmax/1e3/1.4;  % mass from power density ~1.4 kW/kg
mb    = battery_weight;
mtr   = 17;
m0    = 1180;      % i3 curb mass
mcargo= 0;
mv    = m0 + mb + mem + mtr + mcargo;

% i3 rolling/drag
f_r    = 0.4;
lambda = 1.05;
cr     = 0.0174;
cd     = 0.29;
Af     = 2.38;
dw     = 0.6996;

%% Performance calculations
g = 9.81;
rho=1.18;
kR = 0.55;
muP=0.90;

Ttmax = kR*muP*mv*g*dw/2;
Ptmax = Pemmax*0.97;
vb = Ptmax./(Ttmax./(dw/2));
vf = 100/3.6;
ta = (lambda*mv*(vb^2+vf^2))./(2*(Ptmax - (2/3)*mv*g*cr*vf - (1/5)*rho*cd*Af*vf^3));

% Max. speed with 2nd gear => (dw/2)*w_EM_max/g2
vmax = (dw/2)*max(w_EM_max)/g2;

%% Simulation
results = sim(sys_name);
cons_result = results.Consumption(end);

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
    results.t, results.P_EM1/1e3, results.t, (results.P_EM1+results.Losses)/1e3); 
grid;
xlabel('Time [s]');
ylabel('[kW]');
legend('P_{wheel}','P_{trans}', 'P_{EM}', 'P_{BT,output}');

fig = figure;
set(fig,'NumberTitle', 'off')
set(fig,'Name', 'Powers')
plot(results.P_wheel/1e3, results.P_MGB/1e3,'.', ...
    results.P_wheel/1e3, results.P_EM1/1e3,'o', ...
    results.P_wheel/1e3, (results.P_EM1+results.Losses)/1e3,'.');
grid; 
legend('P_{trans}', 'P_{EM}', 'P_{BT,output}');
xlabel('P_{wheel} [kW]'); ylabel('Powers [kW]');

fig = figure;
set(fig,'NumberTitle', 'off') 
set(fig,'Name', 'Battery state of charge')
plot(linspace(min(results.t),max(results.t),length(results.SoC)), results.SoC);
grid;
ylabel('SoC [%]');
xlabel('Time [s]');
ylim([0 100]);

fig = figure;
set(fig,'NumberTitle', 'off')
set(fig,'Name', 'Operation points')
v = 0:0.05:2;
[x,y] = meshgrid(w_EM_row*30/pi, T_EM_col);
[c,h] = contourf(x,y, eta_EM_mapM',v); 
clabel(c,h); 
hold on;
plot(w_EM_max*30/pi, T_EM_max, w_EM_max*30/pi, -T_EM_max, 'r', 'linewidth', 2);
plot(results.w_EM1*30/pi, results.T_EM1, 'o');
xlabel('w_{EM} [rpm]'); ylabel('T_{EM} [Nm]'); 
grid;
