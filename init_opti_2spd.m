%% 4AUB10 - EV 2-Speed Optimization
% Optimize [g1, g2, scale_EM, shift] to minimize energy consumption (kWh/100km)
% Subject to: top speed >= 150 km/h and acceleration time <= 11 s.
% Designs that trigger an overload/overspeed error or violate constraints return Inf.
clear;
global N_sim
N_sim = 1800;  % Number of simulation steps

%% Design Vector Initialization
% 4-element design vector: [g1, g2, scale_EM, shift]
x0 = [10, 6, 4, 28];

%% Optimization Bounds and Objective Function
lb = [5, 1, 1, 10];
ub = [15, 8, 8, 300];
objHandle = @(paramVec) objFun_2spd(paramVec);

%% Run Optimization (Pattern Search)
options = optimoptions('patternsearch', 'Display', 'iter', 'MaxIterations', 40);
[x_best, f_best] = patternsearch(objHandle, x0, [], [], [], [], lb, ub, [], options);

%% Post-Optimization Verification
[f_val_check, overspeed_flag] = objHandle(x_best);
if overspeed_flag || isinf(f_val_check)
    warning('The best design found encountered an overspeed/overload error or did not meet the constraints. Please revise your design constraints.');
else
    disp('A valid design was found that meets all constraints.');
end

%% Display Optimization Results
disp('--- Pattern Search Completed ---');
fprintf('Best design found: [g1 = %.2f, g2 = %.2f, scale_EM = %.2f, shift = %.2f]\n', ...
    x_best(1), x_best(2), x_best(3), x_best(4));
fprintf('Minimum consumption: %.3f [kWh/100km]\n', f_best);
