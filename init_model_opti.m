%% Setting up parameters of different componenets
global g1
global g2
global scale_EM
global Np
global Ns
global s1
global mv
% Gear box 
e_gb    = 0.97;   % internal efficiency [-]
Ploss   = 100;    % stationary losses [W]
wem_min = 1;      % Minimum wheel speed beyond which losses are generated [rad/s]

% Init of the Multispeed Transmission block
path= [sys_name,'/Multispeed gearbox'];
set_param(path,'nr_gears', '2','i_1', 'g1', 'i_2','g2','i_diff','1', 'e_GT', 'e_gb', 'P_GT0', 'Ploss', 'w_wheel_min', 'wem_min' );

% Battery 
init_SoC=90;
path= [sys_name,'/Battery (ECM)'];
set_param(path,"init_SoC",'init_SoC');
set_param(path,"bt_Np",'Np');
set_param(path,"bt_Ns",'Ns');
set_param(path,"SoC_low",'5');
set_param(path,"SoC_high",'95');
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
set_param(path,'P_aux','Paux');
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
clear components;
 
% vehicle mass
mem=Pemmax/1e3/1.4;                 % E-machine mass [kg]
mb     = battery_weight;          % battery mass [kg]
mtr    = 35;                                  % 2spd transmission mass
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

% Init of the Vechile body block
path= [sys_name,'/Vehicle body'];
set_param(path,'m_f','mv','A_f','Af','d_wheel','dw','cw','cd','mu','cr');

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









