clc; clear all; close all;
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath("../Deliverable_3_1"));
rmpath(genpath("../Deliverable_5_1"));
addpath(fullfile('..', 'src'));
addpath(fullfile('..', "Deliverable_3_2"))

%% Set the parameters of the simulations

Ts = 1/10; % Note that we choose a larger Ts here to speed up the simulation
rocket = Rocket(Ts);
H = 0.5;
Tf = 30;

%% With 15° for the roll

% Define the controller
nmpc = NMPC_Control(rocket, H);

% MPC reference with default maximum roll = 15 deg
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
x0 = zeros(12,1);
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, nmpc, ref);

rocket.anim_rate = 10; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
saveas(ph.fig,"Graphs/nmpc_15_roll.svg");

%% With 50° for the roll

% Define the controller 
nmpc = NMPC_Control(rocket, H);

% MPC reference with specified maximum roll = 50 deg
roll_max = deg2rad(50);
ref = @(t_, x_) rocket.MPC_ref(t_, Tf, roll_max);
x0 = zeros(12,1);
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, nmpc, ref);

rocket.anim_rate = 10; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
saveas(ph.fig,"Graphs/nmpc_50_roll.svg");

%% Linear controller with 50° roll

[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);
% Design MPC controller
H = 1.5; % Horizon length in seconds


% Merge four sub−system controllers into one full−system controller
mpc_x = MPC_Control_x(sys_x, Ts, H);
mpc_y = MPC_Control_y(sys_y, Ts, H);
mpc_z = MPC_Control_z(sys_z, Ts, H);
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);

mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

% Setup reference function
roll_max = deg2rad(50);
ref = @(t_, x_) rocket.MPC_ref(t_, Tf, roll_max);
x0 = zeros(12,1);

[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);

% Plot pose
rocket.anim_rate = 10; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title

saveas(ph.fig,"Graphs/Linear_50roll.svg");

%% Absolute cost function - 15° for the roll

H = 0.5;

% Define the controller
nmpc = NMPC_Control_abs_cost(rocket, H);

% MPC reference with default maximum roll = 15 deg
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
x0 = zeros(12,1);
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, nmpc, ref);

rocket.anim_rate = 10; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
saveas(ph.fig,"Graphs/non_linear_nmpc_15_roll.svg");

