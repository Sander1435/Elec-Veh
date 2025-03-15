%% init_opti_2spd.m
% Example of using patternsearch to optimize [g1, g2, scaleEM, shift]
% by calling objFun_2spd(x).

clear; clc;


% Let's define an initial guess for x:
% [g1, g2, scale_EM, shift]
x0 = [10, 6, 4, 28];

% Bounds (example):
%   g1 in [5..15], g2 in [1..8], scale in [1..8], shift in [10..300]
lb = [ 5, 1, 1, 10 ];
ub = [15, 8, 8, 500 ];

objHandle = @(paramVec) objFun_2spd(paramVec);

options = optimoptions('patternsearch','Display','iter','MaxIterations',20);

[x_best, f_best] = patternsearch(objHandle, x0, [], [], [], [], lb, ub, [], options);

disp('--- Pattern Search done ---');
fprintf('Best design found: [g1=%.2f, g2=%.2f, scale=%.2f, s1=%.2f]\n', ...
    x_best(1), x_best(2), x_best(3), x_best(4));
fprintf('Minimum consumption: %.3f [kWh/100km]\n', f_best);
