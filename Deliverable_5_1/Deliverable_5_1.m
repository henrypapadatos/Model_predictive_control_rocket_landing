addpath(fullfile('..', 'src'));
%%
d1 = 0;
d2 = 0;
%  Pavg = (P1 + P2)/2, limited to [20%, 80%]
Pavg = 0;
%Pdiff = P2 − P1, Pdiff might be up to ±20%
Pdiff = 0;
w = [0 0 0];
phi = [0 0 0];
v = [0 0 0];
p = [0 0 0];

Ts = 1/20;
rocket = Rocket(Ts);
u = [d1, d2, Pavg, Pdiff]'; % (Assign appropriately)
[b_F, b_M] = rocket.getForceAndMomentFromThrust(u)
x = [w, phi, v, p]'; % (Assign appropriately)
%derivee de x
x_dot = rocket.f(x, u)

rocket = Rocket(Ts);
Tf = 2.0; % Time to simulate for
x0 = [deg2rad([2 -2 0, -2 2 0]), 0 0 0, 0 0 0]'; % (w, phi, v, p) Initial state
u = [deg2rad([2 0]), 60, 0 ]'; % (d1 d2 Pavg Pdiff) Constant input
[T, X, U] = rocket.simulate_f(x0, Tf, u); % Simulate nonlinear dynamics f
rocket.anim_rate = 1.0;
rocket.vis(T, X, U); 
%%
%%TO DO 5

Ts = 1/20; %Sample time
rocket = Rocket(Ts);
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);
% Design MPC controller
H = 5; % Horizon length in seconds

mpc_z = MPC_Control_z(sys_z, Ts, H);

rocket.mass = 1.783; % Manipulate mass for simulation
z_ref=1;
z = [0;0];

% Get control input
uz = mpc_z.get_u(z,z_ref);

Tf = 10;
z0 = [0;0];
[T, Z_sub, U_sub] = rocket.simulate(sys_z, z0, Tf, @mpc_z.get_u, z_ref);
%[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);


