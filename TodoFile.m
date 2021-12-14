clear all, close all

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

Tf = 0.5;
y0 = [0;0;0;0.5];
[T, Y_sub, U_sub] = rocket.simulate(sys_y, y0, Tf, @mpc_y.get_u, 0);
ph = rocket.plotvis_sub(T, Y_sub, U_sub, sys_y, xs, us);

%% Todo 3.1: confirming good functioning of MPC_Control_z

figure(3)
mpc_z = MPC_Control_z(sys_z, Ts, H);
z = [0.1;0.2];
% Get control input
uz = mpc_z.get_u(z);

Tf = 0.5;
z0 = [0;0.5];
[T, Z_sub, U_sub] = rocket.simulate(sys_z, z0, Tf, @mpc_z.get_u, 0);
ph = rocket.plotvis_sub(T, Z_sub, U_sub, sys_z, xs, us);

%% Todo 3.1: confirming good functioning of MPC_Control_roll

figure(4)
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);
roll = [0.1;0.1];
% Get control input
uroll= mpc_roll.get_u(roll);

Tf = 0.5;
roll0 = [0;0.5];
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
x = [0;0;0.1;0];
x_ref = 1;
% Get control input
ux = mpc_x.get_u(x, x_ref);

Tf = 2;
x0 = [0;0;0;0.5];
[T, X_sub, U_sub] = rocket.simulate(sys_x, x0, Tf, @mpc_x.get_u, x_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us, x_ref);

