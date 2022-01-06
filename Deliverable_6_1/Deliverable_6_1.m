addpath(fullfile('..', 'src'));


%% With 15° for the roll

Ts = 1/10; % Note that we choose a larger Ts here to speed up the simulation
rocket = Rocket(Ts);
H = 0.5;
nmpc = NMPC_Control_roll15(rocket, H);

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

%% With 50° for the roll

Ts = 1/10; % Note that we choose a larger Ts here to speed up the simulation
rocket = Rocket(Ts);
H = 0.5;
nmpc = NMPC_Control_roll50(rocket, H);

% MPC reference with default maximum roll = 15 deg
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
% MPC reference with specified maximum roll = 50 deg
roll_max = deg2rad(50);
ref = @(t_, x_) rocket.MPC_ref(t_, Tf, roll_max);
x0 = zeros(12,1);
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, nmpc, ref);

rocket.anim_rate = 1; % Increase this to make the animation faster
                      % anim rate = 4 is about right for printing in the report
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. sim'; % Set a figure title
