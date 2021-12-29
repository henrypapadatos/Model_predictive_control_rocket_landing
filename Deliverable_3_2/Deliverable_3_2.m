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

mpc_x = MPC_Control_x(sys_x, Ts, H);
x_ref = -5;
x = [0;0;0.1;0.5];
% 
% % Get control input
ux = mpc_x.get_u(x, x_ref);

Tf = 5;
x0 = [0;0;0;0];
[T, X_sub, U_sub] = rocket.simulate(sys_x, x0, Tf, @mpc_x.get_u, x_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us, x_ref);

%% MPC_Control_y

mpc_y = MPC_Control_y(sys_y, Ts, H);
y_ref = 5;
y = [0;0;0.1;0.5];

% Get control input
uy = mpc_y.get_u(y, y_ref);

Tf = 5;
y0 = [0;0;0;0];
[T, Y_sub, U_sub] = rocket.simulate(sys_y, y0, Tf, @mpc_y.get_u, y_ref);
ph = rocket.plotvis_sub(T, Y_sub, U_sub, sys_y, xs, us, y_ref);

%% MPC_Control_z

mpc_z = MPC_Control_z(sys_z, Ts, H);
z_ref= 5;

z = [0;0];

% Get control input
uz = mpc_z.get_u(z,z_ref);

Tf = 5;
z0 = [0;0];
[T, Z_sub, U_sub] = rocket.simulate(sys_z, z0, Tf, @mpc_z.get_u, z_ref);
ph = rocket.plotvis_sub(T, Z_sub, U_sub, sys_z, xs, us);

%% MPC_Control_roll

mpc_roll = MPC_Control_roll(sys_roll, Ts, H);
roll = [0;0.8];
roll_ref = 0.5;
% Get control input
uroll= mpc_roll.get_u(roll, roll_ref);

Tf = 5;
roll0 = [0;0.8];
[T, roll_sub, U_sub] = rocket.simulate(sys_roll, roll0, Tf, @mpc_roll.get_u, roll_ref);
ph = rocket.plotvis_sub(T, roll_sub, U_sub, sys_roll, xs, us);