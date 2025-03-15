% ==============
% TU/e QSS Toolbox
% CVT controller block initialization script
% ==============

%% Global variables
global h                                       % Stepsize [s] from block "Driving Cycle"

%%
%if ~exist('T_EM_col')
%load eta_EM_mapG                    % Efficiency map                [-]

%load w_EM_row                       % Motor speed range             [rad/s]
%load T_EM_col                       % Motor torque range            [Nm]
%load w_EM_max                       % Maximum motor speed           [rad/s]    
%load T_EM_max                       % Maximum motor torque          [Nm]
%end

eta_EM_mapG = evalin('base', 'eta_EM_mapG')';
w_EM_row=evalin('base', 'w_EM_row');
w_EM_max=evalin('base', 'w_EM_max');
T_EM_col=evalin('base', 'T_EM_col');
T_EM_max=evalin('base', 'T_EM_max');

%% Pre-processing quantities
eta_EM_map_flip = flipud(eta_EM_mapG');
%eta_EM_map_flip(eta_EM_map_flip == 10) = 2;

T_EM_pos = T_EM_col(ceil(length(T_EM_col)/2):end);
eta_EM_pos = flip(transpose(eta_EM_mapG(1:ceil(length(T_EM_col)/2),:)));
P_EM = w_EM_row .* T_EM_col(ceil(length(T_EM_col)/2):end);

eta_EM_neg = eta_EM_map_flip(:,1:floor(length(T_EM_col)/2));
T_EM_neg = T_EM_col(1:floor(length(T_EM_col)/2));

%% Compute optmimal operation line
w = 0:10:w_EM_max(end);

% positive
T_p = interp1(w_EM_max,T_EM_max,w);
for i=1:length(w)
    for j = 1:length(T_p)
        eta_local(j) = interp2(w_EM_row,T_EM_pos,eta_EM_pos',w(i),T_p(j));
    end
    [eta_opt(i),ind_T(i)] = max(eta_local);
end
P_opt_pos = w.*T_p(ind_T);

ind = find(w == 300);
w1 = w(1:ind);
w2 = w(ind+1:end);

P_opt_pos1 = P_opt_pos(1:ind);
[P_opt_pos2, i] = sort(P_opt_pos(ind+1:end));
w2S_p = w2(i);

% negative
T_n = interp1(w_EM_max,-T_EM_max,w);
for i=1:length(w)
    for j = 1:length(T_n)
        eta_local(j) = interp2(w_EM_row,T_EM_neg,eta_EM_neg',w(i),T_n(j));
    end
    [eta_opt_n(i),ind_T_n(i)] = min(eta_local);
end
P_opt_neg = w.*T_n(ind_T_n);

ind = find(w == 300);
w1_n = flip(w1);
P_opt_neg1 = flip(P_opt_neg(1:ind));
[P_opt_neg2, i] = sort(P_opt_neg(ind+1:end));
w2S_n = w2(i);

