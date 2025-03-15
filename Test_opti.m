clear all;
%% Global variables
global N_sim
global g1 scale_EM
global var
global ta mv x
global x_result
global y_result
global z_result
%% add QSS toolbox to the path
p = genpath('QSS_TB');
addpath(p);
sys_name = 'EV_1spd_4AUB10_2023b';
load_system(sys_name);
%% Initialization of the optimization routine
x_result = []; y_result = []; z_result = [];
x_guess  = [9;4                               % gear ratio [-], motor scale [-]
x = [x_guess(1) x_guess(2) 24 144 2 2 2 1];  % set of input design parameters, x
g1              = x(1);       % gear ratio [-]
scale_EM        = x(2);       % EM scale [-]
Np              = x(3);       % Number of cells connected in parallel [-]
Ns              = x(4);       % Number of cells connected in series [-]
motortype       = x(5);       % Motortype 1 is a PMSM and Motortype 2 is a IM [-] 
motornumberinit = x(6);       % This parameter changes between 2 powers of an engine
BatterySort     = x(7);       % Decides between High Power(2) or Energy cells(1) [-]
CellSort        = x(8);       % Decides the exact cells used for HE from 1 to 2[-] and for HP 1 to 4[-]

%% Optimization routine
options = optimset('PlotFcns',@optimplotfval);
init_model_opti;
disp('Optimization started...')
disp(' ')
[x_opt, V_result, overspeed] = fminsearch(@Opti_EV_PT_1spd, x_guess, options); % Nelder-Mead simplex method
disp(' ')
disp('...optimization finished!')

%% Plots
figure;
plot(x_result); legend('g1 [-] ','scale_{EM} [-]');
xlabel('Iterations')
ylabel('Optimization values [-]')
title('Variations of the Optimization values during the optimization')

figure;
plot(y_result);
xlabel('Iterations')
ylabel('Vehicle mass (output of model) [kg]')
title('Variations of the Optimization values during the optimization')

figure;
plot(z_result);
xlabel('Iterations')
ylabel('Acceleration time (output of model) [s]')
title('Variations of the Optimization values during the optimization')
