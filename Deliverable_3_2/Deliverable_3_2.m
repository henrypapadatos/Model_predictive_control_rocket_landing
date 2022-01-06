clc; clear all; close all;
addpath(fullfile('..', 'src'));

%% Linearize and split the system

Ts = 1/20; % Sample time
H = 1.5; % Horizon length in seconds
rocket = Rocket(Ts);
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

%% MPC_Control_x

% Initialize the controller
mpc_x = MPC_Control_x(sys_x, Ts, H);

x_ref = -5; % Tracks a reference of -5m
Tf = 8; % Set the duration of the simulation
x0 = [0;0;0;0]; % Starts at origin

[T, X_sub, U_sub] = rocket.simulate(sys_x, x0, Tf, @mpc_x.get_u, x_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us, x_ref);

%% MPC_Control_y

% Initialize the controller
mpc_y = MPC_Control_y(sys_y, Ts, H);

y_ref = -5; % Tracks a reference of -5m
Tf = 8; % Set the duration of the simulation
y0 = [0;0;0;0]; % Starts at origin

[T, Y_sub, U_sub] = rocket.simulate(sys_y, y0, Tf, @mpc_y.get_u, y_ref);
ph = rocket.plotvis_sub(T, Y_sub, U_sub, sys_y, xs, us, y_ref);

%% MPC_Control_z

% Initialize the controller
mpc_z = MPC_Control_z(sys_z, Ts, H);

z_ref = -5; % Tracks a reference of -5m
Tf = 8; % Set the duration of the simulation
z0 = [0;0]; % Starts at origin

[T, Z_sub, U_sub] = rocket.simulate(sys_z, z0, Tf, @mpc_z.get_u, z_ref);
ph = rocket.plotvis_sub(T, Z_sub, U_sub, sys_z, xs, us);

%% MPC_Control_roll

% Initialize the controller
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);

roll_ref = pi/4; % Tracks a reference of 45°
Tf = 8; % Set the duration of the simulation
roll0 = [0;0]; % Starts at origin

[T, roll_sub, U_sub] = rocket.simulate(sys_roll, roll0, Tf, @mpc_roll.get_u, roll_ref);
ph = rocket.plotvis_sub(T, roll_sub, U_sub, sys_roll, xs, us);