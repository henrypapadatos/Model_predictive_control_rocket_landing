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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

f_discrete = @(x,u) RK4(x,u,rocket.Ts,rocket.f);

opti.minimize(sum((X_sym(10,:)-ref_sym(1)).^2) +...
    sum((X_sym(11,:)-ref_sym(2)).^2)+...
    sum((X_sym(12,:)-ref_sym(3)).^2)+...
    sum((X_sym(6,:)-ref_sym(4))).^2);


for k=1:N-1 % loop over control intervals
    rocket.f(X_sym(:,k), U_sym(:,k));
   opti.subject_to(X_sym(:,k+1) == f_discrete(X_sym(:,k), U_sym(:,k)));
%     opti.subject_to(X_sym(:,k+1) == RK4(X_sym(:,k), U_sym(:,k),rocket.Ts,rocket.f));

end

%input constraints
opti.subject_to(-0.26<=U_sym(1)<0.26);
opti.subject_to(-0.26<=U_sym(2)<0.26);
opti.subject_to(50<=U_sym(3)<=80);
opti.subject_to(-20<=U_sym(4)<=20);

%state constraints
opti.subject_to(-deg2rad(85)<=X_sym(5)<=deg2rad(85));

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
