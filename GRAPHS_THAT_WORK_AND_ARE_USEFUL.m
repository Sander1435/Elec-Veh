%%
% --- Section 1: Motor Speed vs. Time with Gear Coloring ---
% Use results.w_EM1 and the defined shift point s1 (in rad/s)
% Assume results.t is the time vector

gear1_idx = results.w_EM1 < s1;  % operating in gear 1
gear2_idx = results.w_EM1 >= s1; % operating in gear 2

figure;
hold on;
plot(results.t(gear1_idx), results.w_EM1(gear1_idx), 'bo', 'DisplayName', 'Gear 1');
plot(results.t(gear2_idx), results.w_EM1(gear2_idx), 'rs', 'DisplayName', 'Gear 2');
xlabel('Time [s]');
ylabel('Motor Speed [rad/s]');
title('Motor Speed vs. Time (Gear 1 vs. Gear 2)');
legend('Location','best');
grid on;
hold off;

%% --- Section 2: Time and Energy Analysis per Gear ---
% Use results.t, results.P_EM1, and the shift point s1 from the simulation

dt = [diff(results.t); 0];  % approximate time step at each sample

% Define gear masks based on the motor speed (w_EM1) compared to s1.
gear1_mask = results.w_EM1 < s1;
gear2_mask = results.w_EM1 >= s1;

% Total time spent in each gear
timeGear1 = sum(dt(gear1_mask));
timeGear2 = sum(dt(gear2_mask));

% Compute energy (in Joules) from motor power (P_EM1, in Watts)
energyGear1 = sum(results.P_EM1(gear1_mask) .* dt(gear1_mask));
energyGear2 = sum(results.P_EM1(gear2_mask) .* dt(gear2_mask));

fprintf('Time in Gear 1: %.2f s\n', timeGear1);
fprintf('Time in Gear 2: %.2f s\n', timeGear2);
fprintf('Energy in Gear 1: %.2f J\n', energyGear1);
fprintf('Energy in Gear 2: %.2f J\n', energyGear2);



%% --- Section 3: Detect Shift Events ---
% Find indices where motor speed crosses the shift threshold s1.
% We look for a change from below s1 to above s1 or vice versa.

shiftEvents = find(diff(results.w_EM1 < s1) ~= 0);
shiftTimes = results.t(shiftEvents);

% Plot motor speed with shift event markers
figure;
plot(results.t, results.w_EM1, 'b-', 'LineWidth', 1.5);
hold on;
plot(shiftTimes, results.w_EM1(shiftEvents), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Motor Speed [rad/s]');
title('Motor Speed with Detected Shift Events');
legend('Motor Speed','Shift Events','Location','best');
grid on;
hold off;

%% --- Section 4: Motor Torque vs. Motor Speed Colored by Gear ---
% Use results.T_EM1 for motor torque and results.w_EM1 for speed

gear1_mask = results.w_EM1 < s1;
gear2_mask = results.w_EM1 >= s1;

% Convert motor speed from rad/s to rpm for plotting clarity
motorSpeed_rpm = results.w_EM1 * (30/pi);

figure;
hold on;
scatter(motorSpeed_rpm(gear1_mask), results.T_EM1(gear1_mask), 10, 'b', 'filled', 'DisplayName', 'Gear 1');
scatter(motorSpeed_rpm(gear2_mask), results.T_EM1(gear2_mask), 10, 'r', 'filled', 'DisplayName', 'Gear 2');
xlabel('Motor Speed [rpm]');
ylabel('Motor Torque [Nm]');
title('Motor Torque vs. Speed (Gear 1 vs. Gear 2)');
legend('Location','best');
grid on;
hold off;

%% --- Extra Technical Section 1: Wheel Power vs. Road Load Power ---
% Compute vehicle speed at the wheel (m/s):
v_vehicle = (dw/2) * results.w_wheel;  % v = (wheel radius)*w_wheel

% Compute theoretical road load forces on a flat road:
% Rolling resistance force: F_roll = cr * mv * g
F_roll = cr * mv * g;  % constant on flat road

% Aerodynamic drag force: F_aero = 0.5 * rho * cd * Af * v^2
F_aero = 0.5 * rho * cd * Af .* (v_vehicle.^2);

% Total theoretical tractive force:
F_total_theoretical = F_roll + F_aero;

% Theoretical power demand: P_theoretical = F_total_theoretical * v_vehicle (W)
P_theoretical = F_total_theoretical .* v_vehicle;

% Plot simulated wheel power (P_wheel) and theoretical power vs. time:
figure;
plot(results.t, results.P_wheel/1e3, 'b-', 'LineWidth', 2); hold on;
plot(results.t, P_theoretical/1e3, 'r--', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Power [kW]');
title('Simulated Wheel Power vs. Theoretical Road Load Power');
legend('Simulated P_{wheel}', 'Theoretical P_{road}');
grid on;
hold off;
