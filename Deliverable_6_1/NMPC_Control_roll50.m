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

% Define cost matrices
Q = eye(nx);
Q(1,1) = 600; % Minimize angular velocity about x
Q(2,2) = 1100; % Minimize angular velocity about y
Q(3,3) = 50; % Minimize angular velocity about z
Q(4,4) = 1; % Minimize alpha angle
Q(5,5) = 1; % Minimize beta angle
Q(6,6) = 400; % Track roll reference 
Q(7,7) = 1; % Minimize velocity about x
Q(8,8) = 1; % Minimize velocity about y
Q(9,9) = 1; % Minimize velocity about z
Q(10,10) = 1000; % Track x reference
Q(11,11) = 1000; % Track y reference
Q(12,12) = 1000; % Track z reference

% Define cost matrix for the input
R = [1 0 0 0; % Minimize delta 1
     0 1 0 0; % Minimize delta 2
     0 0 5 0; % Minimize P_avg
     0 0 0 0.5]; % Minimize P_diff

% Linearize the system to compute terminal cost

X0 = SX.sym('X0',nx,1); % declare a symbolic variable X0 of size nxx1
U0 = SX.sym('U0',nu,1); % declare a symbolic variable U0 of size nux1

x_next = rocket.f(X0,U0); % Calculate the next step symbolically

A_jac = jacobian(x_next,X0); % returns jacobian for an expression (x_next) w.r.t X0
B_jac = jacobian(x_next,U0); % returns jacobian for an expression (x_dot) w.r.t U0

% convert A_algorithmic, B_algorithmic expressions to a callable functions
% that will return jacobian w.r.t x and w.r.t. u
A_algorithmic = casadi.Function('A_algorithmic',{X0,U0},{A_jac});
B_algorithmic = casadi.Function('A_algorithmic',{X0,U0},{B_jac});

% Equilibrum point around which the system will be linearized
[xs, us] = rocket.trim();

% Discretize the system
A = eye(nx) + rocket.Ts*full(A_algorithmic(xs,us));
B = rocket.Ts*full(B_algorithmic(xs,us));

% Compute the terminal cost matrix
[~,P,~] = dlqr(A,B,Q,R);

% Define the target tracked by the rocket
x_target = [0 0 0 0 0 ref_sym(4) 0 0 0 ref_sym(1:3)']';
u_target = us;

% Define the matrices for the augmented system
Q_aug = blkdiag(kron(eye(N-1),Q),P);
R_aug = kron(eye(N-1),R);

% Discretize f
f_discrete = @(x,u) RK4(x,u,rocket.Ts,@rocket.f);

% Define the optimization problem
opti.minimize((reshape(X_sym - x_target,nx*N,1)')*Q_aug*reshape(X_sym - x_target,nx*N,1)+...
              (reshape(U_sym - u_target,nu*(N-1),1)')*R_aug*reshape(U_sym - u_target,nu*(N-1),1)...
              );

% Define the constraints 

% Dynamic constraints 
for k=1:N-1 % loop over control intervals
    opti.subject_to(X_sym(:,k+1) == f_discrete(X_sym(:,k), U_sym(:,k)));
end

% Input constraints
opti.subject_to(-0.26 <= U_sym(1) <= 0.26);
opti.subject_to(-0.26 <= U_sym(2) <= 0.26);
opti.subject_to(50 <= U_sym(3) <= 80);
opti.subject_to(-20 <= U_sym(4) <= 20);

% State constraints
opti.subject_to(-deg2rad(85) <= X_sym(5) <= deg2rad(85));

% Initial state constraints
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
