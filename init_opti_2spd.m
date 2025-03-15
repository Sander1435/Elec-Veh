%% init_opti_2spd.m
% Example of how to call the objective function objFun_2spd.m 
% with a built-in optimizer

clear; clc;

% Possibly define a global "N_sim" if your model uses that:
global N_sim
N_sim = 3000;  % or whatever your drive cycle final time is in steps

% Let's guess an initial design:
%   Suppose x = [g1, g2, scale_EM, Np, Ns, motortype, motornumberinit, BatterySort, CellSort, s1]
x0 = [12, 6, 4, 16, 96, 1, 1, 1, 1, 300];

% Lower/upper bounds:
lb = [ 5,  3, 1,  1,  1, 1, 1, 1, 1, 100]; 
ub = [15, 10,10, 50,200, 2, 2, 2, 4,1000];

% Define the objective function handle:
objHandle = @(xx) objFun_2spd(xx);

% Pick an optimizer, e.g. "patternsearch" or "ga"
opts = optimoptions('patternsearch','Display','iter','MaxIterations',30);

% Run the optimization
[x_best, fval_best] = patternsearch(objHandle, x0, [], [], [], [], lb, ub, [], opts);

% Show final results
fprintf('Best design found:\n');
disp(x_best);
fprintf('Minimum consumption: %g [kWh/100km]\n', fval_best);
