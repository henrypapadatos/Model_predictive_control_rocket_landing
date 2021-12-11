clc;
close all;
clear all;

addpath('./src')

Ts = 1/20;
rocket = Rocket(Ts);

[xs, us] = rocket.trim() % Compute steady−state for which 0 = f(xs,us)
sys = rocket.linearize(xs, us) % Linearize the nonlinear model about trim point