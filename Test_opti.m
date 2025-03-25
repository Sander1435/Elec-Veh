clear all;
%% Global variables
global N_sim
global g1
global g2
global scale_EM
global Np
global Ns
global s1
global var
global ta mv x
global x_result
global y_result
global z_result
%% add QSS toolbox to the path
p = genpath('F:\Github\Elec-Veh');
addpath(p);
sys_name = 'EV_2spd_4AUB10_2023b';
load_system(sys_name);
%% Initialization of the optimization routine

x_result = []; y_result = []; z_result = [];
x_guess  = [10;4;2.214;24;144;100]; 
x = [x_guess(1) x_guess(2) x_guess(3) x_guess(4) x_guess(5) 2 2 1 2 x_guess(6)];  % set of input design parameters, x
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

%% Optimization routine
options = optimset('PlotFcns',@optimplotfval);
init_model_opti;
disp('Optimization started...')
disp(' ')
[x_opt, V_result, overspeed] = fminsearch(@Opti_EV_PT_2spd, x_guess, options); % Nelder-Mead simplex method
disp(' ')
disp('...optimization finished!')

%% Plots
figure;
plot(x_result); legend('g1 [-] ','g2 [-]','scale_EM [-]', 'Np [-]', 'Ns [-]', 's1 [-]');
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
