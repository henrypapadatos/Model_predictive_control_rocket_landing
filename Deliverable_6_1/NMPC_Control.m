function opti_eval = NMPC_Control(rocket, H)

import casadi.*
opti = casadi.Opti(); % Optimization problem

N = ceil(H/rocket.Ts); % MPC horizon
nx = 12; % Number of states
nu = 4;  % Number of inputs

% Decision variables (symbolic)
X_sym = opti.variable(nx, N); % state trajectory
U_sym = opti.variable(nu, N-1);   % control trajectory)

% Parameters (symbolic)
x0_sym  = opti.parameter(nx, 1);  % initial state
ref_sym = opti.parameter(4, 1);   % target position

% Slack variables 
%epsilon_beta_1 = opti.variable(1,N);
%epsilon_beta_2 = opti.variable(1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

f_discrete = @(x,u) RK4(x,u,rocket.Ts,@rocket.f);

opti.minimize(...
    5*sum((X_sym(10,:)-ref_sym(1)).^2) +... % Track x reference
    5*sum((X_sym(11,:)-ref_sym(2)).^2) +... % Track y reference
    5*sum((X_sym(12,:)-ref_sym(3)).^2) +... % Track z reference
    40*sum((X_sym(6,:)-ref_sym(4)).^2) +... % Track roll reference
    30*sum(X_sym(1,:).^2) +... % Minimize angular velocity about x
    60*sum(X_sym(2,:).^2) +... % Minimize angular velocity about y
    0*sum(X_sym(3,:).^2) +... % Minimize angular velocity about z
    0*sum(X_sym(7,:).^2) +... % Minimize velocity along x
    0*sum(X_sym(8,:).^2) +... % Minimize velocity along y
    0*sum(X_sym(9,:).^2) +... % Minimize velocity along z
    0.02*sum((U_sym(3,:) - 60).^2) + ... % Minimize the energy used 
    0.25*sum(U_sym(4,:).^2)... % Minimize P_diff
    );

%    epsilon_beta_1*epsilon_beta_1' + sum(epsilon_beta_1) +... % Slack constraints
%    epsilon_beta_2*epsilon_beta_2' + sum(epsilon_beta_2));


for k=1:N-1 % loop over control intervals
    opti.subject_to(X_sym(:,k+1) == f_discrete(X_sym(:,k), U_sym(:,k)));
end

%input constraints
opti.subject_to(-0.26 <= U_sym(1) <= 0.26);
opti.subject_to(-0.26 <= U_sym(2) <= 0.26);
opti.subject_to(50 <= U_sym(3) <= 80);
opti.subject_to(-20 <= U_sym(4) <= 20);

%state constraints
opti.subject_to(-deg2rad(85) <= X_sym(5) <= deg2rad(85));

%initialize state
opti.subject_to(X_sym(:,1) == x0_sym);

% YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Setup solver ------
ops = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
opti.solver('ipopt', ops);

% Create function to solve and evaluate opti
opti_eval = @(x0_, ref_) solve(x0_, ref_, opti, x0_sym, ref_sym, U_sym);
end

function u = solve(x0, ref, opti, x0_sym, ref_sym, U_sym)

% ---- Set the initial state and reference ----
opti.set_value(x0_sym, x0);
opti.set_value(ref_sym, ref);

% ---- Solve the optimization problem ----
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');

u = opti.value(U_sym(:,1));

% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end
