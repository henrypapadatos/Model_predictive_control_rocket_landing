clc;

addpath('./src')

Ts = 1/20;
rocket = Rocket(Ts);

Tf = 2.0; % Time to simulate for
% x0 = [deg2rad([2 -2 0, -2 2 0]), 0 0 0, 0 0 0]'; % (w, phi, v, p) Initial state
% u = [deg2rad([2 0]), 60, 0 ]'; % (d1 d2 Pavg Pdiff) Constant input

x0 = [deg2rad([0 0 0, 0 0 0]), 0 0 0, 0 0 0]'; % (w, phi, v, p) Initial state
u = [deg2rad([0 0]), 60, 60 ]'; % (d1 d2 Pavg Pdiff) Constant input

[T, X, U] = rocket.simulate_f(x0, Tf, u); % Simulate nonlinear dynamics f
rocket.anim_rate = 0.5;
rocket.vis(T, X, U); % Trajectory visualization at 1.0x real−time
