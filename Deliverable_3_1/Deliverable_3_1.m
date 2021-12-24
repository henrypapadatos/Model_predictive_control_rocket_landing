clc; clear all; close all;
addpath(fullfile('..', 'src'));
addpath("Graphs")

%% Linearize and split the system

Ts = 1/20;
rocket = Rocket(Ts);
[xs, us] = rocket.trim(); % Compute steady−state for which 0 = f(xs,us)
sys = rocket.linearize(xs, us); % Linearize the nonlinear model about trim point
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

%% 1st system (sys_x [wy beta vx x]): Initialize the controller

% Choose the horizon length
H = 1.1;
% Initialize the controller
mpc_x = MPC_Control_x(sys_x, Ts, H);

%% 1st system (sys_x [wy beta vx x]): Perform simulation

Tf = 8; % Set the duration of the simulation
x0 = [0;0;0;5]; % Starting stationnary at five meters from the origin
[T, X_sub, U_sub] = rocket.simulate(sys_x, x0, Tf, @mpc_x.get_u, 0);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us);
saveas(ph.fig,"Graphs/sys_x_results.svg");

%% 2nd system (sys_y [wx alpha vy y]): Initialize the controller

% Choose the horizon length
H = 1.2;
% Initialize the controller
mpc_y = MPC_Control_y(sys_y, Ts, H);

%% 2nd system (sys_y [wx alpha vy y]): Perform simulation

Tf = 8; % Set the duration of the simulation
x0 = [0;0;0;5]; % Starting stationnary at five meters from the origin
[T, X_sub, U_sub] = rocket.simulate(sys_y, x0, Tf, @mpc_y.get_u, 0);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_y, xs, us);
saveas(ph.fig,"Graphs/sys_y_results.svg");

%% 3rd system (sys_z [vz z]): Initialize the controller

% Choose the horizon length
H = 1.5;
% Initialize the controller
mpc_z = MPC_Control_z(sys_z, Ts, H);

%% 3rd system (sys_z [vz z]): Perform simulation

Tf = 8; % Set the duration of the simulation
x0 = [0;5]; % Starting stationnary at five meters from the origin
[T, X_sub, U_sub] = rocket.simulate(sys_z, x0, Tf, @mpc_z.get_u, 0);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_z, xs, us);
saveas(ph.fig,"Graphs/sys_z_results.svg");

%% 4th system (sys_roll [wz gamma]): Initialize the controller

% Choose the horizon length
H = 0.5;
% Initialize the controller
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);

%% 4th system (sys_roll [wz gamma]): Perform simulation

Tf = 8; % Set the duration of the simulation
x0 = [0;deg2rad(45)]; % Starting stationnary at five meters from the origin
[T, X_sub, U_sub] = rocket.simulate(sys_roll, x0, Tf, @mpc_roll.get_u, 0);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_roll, xs, us);
saveas(ph.fig,"Graphs/sys_roll_results.svg");
