%% Test procedure  - gradient from standstill
%
% Required inputs
%  1. Workspace with loaded model parameters
%  2. Model name as string in variable 'sys_name'
%  Important: driving cycle block should be named 'Driving cycle'
%  Important: model should run without other errors, overloads and etc.
%
% Output
%  1. Vehicle top speed in [rad] in variable 'gradStand_test'
%  Important: result is viable only if overload message is displayed in the
%  command window

%% Load system
load_system('EV_2spd_4AUB10_2023b');

%% Save current cycle 
old_cycle = get_param('EV_2spd_4AUB10_2023b'+ "Driving cycle",'3');

%% Change driving cycle
set_param('EV_2spd_4AUB10_2023b'+ "/Driving cycle",'cyclenr','3');

%% Run model
results_gradSt = sim(sys_name,'SrcWorkspace','current');

%% Determine top speed
gradStand_test = results_gradSt.slope(end);

%% Load old driving cycle 
set_param(sys_name+ "/Driving cycle",'cyclenr',old_cycle); 

%% Clear data
clearvars results_gradSt;