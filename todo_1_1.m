clc;

addpath('./src')

Ts = 1/20;
rocket = Rocket(Ts);

d1=0;
d2=0;
Pavg=0;
Pdiff=0;
w=[0,0,0];
phi=[0,0,0];
v=[0,0,10];
p=[0,0,0];

u = [d1, d2, Pavg, Pdiff]'; % (Assign appropriately)
[b_F, b_m] = rocket.getForceAndMomentFromThrust(u)
x = [w, phi, v, p]'; % (Assign appropriately)
x_dot = rocket.f(x, u)
