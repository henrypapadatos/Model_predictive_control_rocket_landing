clear all, close all
clc

addpath('./src')



%% Gather all the Todo steps of the project

%% Todo 1.1

d1 = 0;
d2 = 0;
Pavg = 0;
Pdiff = 0;
w = [0,0,0];
phi = [0,0,0];
v = [0,0,0];
p = [0,0,0];

Ts = 1/100;
rocket = Rocket(Ts);
u = [d1, d2, Pavg, Pdiff]'; % (Assign appropriately)
[b_F, b_M] = rocket.getForceAndMomentFromThrust(u);
x = [w, phi, v, p]'; % (Assign appropriately)
x_dot = rocket.f(x, u);

%% Todo 1.2

rocket = Rocket(Ts);
Tf = 2.0; % Time to simulate for
% x0 Given in the homework
x0 = [deg2rad([2 -2 0, -2 2 0]), 0 0 0, 0 0 0]'; % (w, phi, v, p) Initial state
% New x0 more convenient at first
x0_easy = [deg2rad([0 0 0, 0 0 0]), 0 0 0, 0 0 0]'; % (w, phi, v, p) Initial state with 0 values
u_ascent = [deg2rad([0 0]), 60, 0 ]'; % Ascending vertically (Translation over z+)
u_descent = [deg2rad([180 0]), 0, 0 ]'; % Descending vertically (Translation over z-)
u_translate_x = [deg2rad([0 -0.1]), 60, 0 ]'; % Translation over x+
u_translate_y = [deg2rad([-0.1 0]), 60, 0 ]'; % Translation over y+
u_rotate_z = [deg2rad([0 0]), 60,  30]'; % Rotation about z
% u_rotate_y = [deg2rad([0 5]), 60,  30]'; % Rotation about y
% u_rotate_x = [deg2rad([0 0]), 60,  30]'; % Rotation about x

% Let see what gives the simulation
[T, X, U] = rocket.simulate_f(x0, Tf, u_ascent); % Simulate nonlinear dynamics f
rocket.anim_rate = 1.0;
rocket.vis(T, X, U); % Trajectory visualization at 1.0x real−time

%% Todo 2.1

rocket = Rocket(Ts);
[xs, us] = rocket.trim(); % Compute steady−state for which 0 = f(xs,us)
sys = rocket.linearize(xs, us); % Linearize the nonlinear model about trim point

% Here is the description for trim function

% Since we know here that a trim point f(xs, us) = 0 exists,
% we can simply solve an box-constrained optimization problem
% with Matlab's fmincon solver:
%                   min f(xs,us).^2
%
% Given that we have nonlinear dynamics, we will find the local
% minimum next to the initial guess.
%
% Decision variable y = [x; u]
            
% For the linearize function it linearizes around the equilibrium point
% using Taylor expansion as follows:
% f(x,u) = f(xs,us) + J_f|x(xs)*(x-xs) + J_f|u(us)*(u-us)

% Let's have a look at the matrices 
A = sys.A; 
B = sys.B;
C = sys.C;
D = sys.D;

% Visualize block matrices now
new_order = [11,8,4,1,10,7,5,2,12,9,6,3];
A_reorder = A(:,new_order);
A_reorder = A_reorder(new_order,:);
B_reorder = B(new_order,:);

%% Todo 2.2

[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

%% Todo 3.1
rmpath("./Deliverable_3_2")
addpath('./Deliverable_3_1')
Ts = 1/20; % Sample time
rocket = Rocket(Ts);
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);
% Design MPC controller
H = 2; % Horizon length in seconds

%% Todo 3.1: confirming good functioning of MPC_Control_x

figure(1)
mpc_x = MPC_Control_x(sys_x, Ts, H);
x = [0;0;0.1;0.5];
% Get control input
ux = mpc_x.get_u(x);

Tf = 4;
x0 = [0;0;0;0.5];
[T, X_sub, U_sub] = rocket.simulate(sys_x, x0, Tf, @mpc_x.get_u, 0);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us);

%% Todo 3.1: confirming good functioning of MPC_Control_y

figure(2)
mpc_y = MPC_Control_y(sys_y, Ts, H);
y = [0;0;0.1;0.5];
% Get control input
uy = mpc_y.get_u(y);

Tf = 5;
y0 = [0;0;0;0.5];
[T, Y_sub, U_sub] = rocket.simulate(sys_y, y0, Tf, @mpc_y.get_u, 0);
ph = rocket.plotvis_sub(T, Y_sub, U_sub, sys_y, xs, us);

%% Todo 3.1: confirming good functioning of MPC_Control_z

figure(3)
mpc_z = MPC_Control_z(sys_z, Ts, H);
z = [0.1;0.2];
% Get control input
uz = mpc_z.get_u(z);

Tf = 4;
z0 = [0;0.5];
[T, Z_sub, U_sub] = rocket.simulate(sys_z, z0, Tf, @mpc_z.get_u, 0);
% ph = rocket.plotvis_sub(T, Z_sub, U_sub, sys_z, xs, us);

%% Todo 3.1: confirming good functioning of MPC_Control_roll

figure(4)
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);
roll = [0;0.8];
% Get control input
uroll= mpc_roll.get_u(roll);

Tf = 5;
roll0 = [0;0.8];
[T, roll_sub, U_sub] = rocket.simulate(sys_roll, roll0, Tf, @mpc_roll.get_u, 0);
ph = rocket.plotvis_sub(T, roll_sub, U_sub, sys_roll, xs, us);

%% Todo 3.2

rmpath("Deliverable_3_1")
addpath("Deliverable_3_2")

Ts = 1/20; % Sample time
rocket = Rocket(Ts);
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);
% Design MPC controller
H = 5; % Horizon length in seconds

%% Todo 3.2: confirming good functioning of MPC_Control_x

mpc_x = MPC_Control_x(sys_x, Ts, H);
x_ref = -5;
% x = [0;0;0.1;0.5];
% 
% % Get control input
% ux = mpc_x.get_u(x, x_ref);

Tf = 5;
x0 = [0;0;0;0];
[T, X_sub, U_sub] = rocket.simulate(sys_x, x0, Tf, @mpc_x.get_u, x_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us, x_ref);
% plot(T,X_sub(4,:));
%% Todo 3.2: confirming good functioning of MPC_Control_y

mpc_y = MPC_Control_y(sys_y, Ts, H);
y_ref = 5;
y = [0;0;0.1;0.5];

% Get control input
uy = mpc_y.get_u(y, y_ref);

Tf = 5;
y0 = [0;0;0;0];
[T, Y_sub, U_sub] = rocket.simulate(sys_y, y0, Tf, @mpc_y.get_u, y_ref);
ph = rocket.plotvis_sub(T, Y_sub, U_sub, sys_y, xs, us, y_ref);

%% Todo 3.2: confirming good functioning of MPC_Control_z

mpc_z = MPC_Control_z(sys_z, Ts, H);
z_ref= 5;

z = [0;0];

% Get control input
uz = mpc_z.get_u(z,z_ref);

Tf = 5;
z0 = [0;0];
[T, Z_sub, U_sub] = rocket.simulate(sys_z, z0, Tf, @mpc_z.get_u, z_ref);
ph = rocket.plotvis_sub(T, Z_sub, U_sub, sys_z, xs, us);
% figure();
% plot(T,Z_sub(2,:));
% title('z subsystem');

%% Todo 3.2: confirming good functioning of MPC_Control_roll

mpc_roll = MPC_Control_roll(sys_roll, Ts, H);
roll = [0;0.8];
roll_ref = 0.5;
% Get control input
uroll= mpc_roll.get_u(roll, roll_ref);

Tf = 5;
roll0 = [0;0.8];
[T, roll_sub, U_sub] = rocket.simulate(sys_roll, roll0, Tf, @mpc_roll.get_u, roll_ref);
ph = rocket.plotvis_sub(T, roll_sub, U_sub, sys_roll, xs, us);
% figure();
% plot(T,roll_sub(2,:));
% title('roll subsystem');
%% Todo 4
rmpath("Deliverable_3_1")
addpath("Deliverable_3_2")

Ts = 1/20; % Sample time
rocket = Rocket(Ts);
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);
% Design MPC controller
H = 3; % Horizon length in seconds


% Merge four sub−system controllers into one full−system controller
mpc_x = MPC_Control_x(sys_x, Ts, H);
mpc_y = MPC_Control_y(sys_y, Ts, H);
mpc_z = MPC_Control_z(sys_z, Ts, H);
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);

mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

% Setup reference function
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
x0 = zeros(12,1);

[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);
% Plot pose
rocket.anim_rate = 4; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
%% TO DO 5.1
rmpath("Deliverable_3_2")
addpath("Deliverable_5_1")
%% TODO: This file should produce all the plots for the deliverable


Ts = 1/20; % Sample time
rocket = Rocket(Ts);
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
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
x0 = zeros(12,1);
rocket.mass = 1.783;
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);
% Plot pose
rocket.anim_rate = 10; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
%%
%%TO DO 5
rmpath("Deliverable_3_2")
addpath("Deliverable_5_1")

Ts = 1/20; % Sample time
rocket = Rocket(Ts);
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
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
%ref = @(t_, x_) [0;0;t_;0];
x0 = zeros(12,1);
rocket.mass = 1.783;
%[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);
[T, X, U, Ref, Z_hat] = rocket.simulate_f_est_z(x0, Tf, mpc, ref, mpc_z, sys_z);
rocket.anim_rate = 10; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
%% part 6
addpath('Deliverable_6_1');
Ts = 1/10; % Note that we choose a larger Ts here to speed up the simulation
rocket = Rocket(Ts);
H = 2;
nmpc = NMPC_Control(rocket, H);
% MPC reference with default maximum roll = 15 deg
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
% MPC reference with specified maximum roll = 50 deg
roll_max = deg2rad(15);
ref = @(t_, x_) rocket.MPC_ref(t_, Tf, roll_max);
x0 = zeros(12,1);
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, nmpc, ref);

rocket.anim_rate = 1; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
