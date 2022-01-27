classdef MPC_Control_z < MPC_Control

    properties
        A_bar, B_bar, C_bar % Augmented system for disturbance rejection
        L                   % Estimator gain for disturbance rejection
    end
    
    methods
        function mpc = MPC_Control_z(sys, Ts, H)
            mpc = mpc@MPC_Control(sys, Ts, H);
            
            [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
        end
        
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   d_est        - disturbance estimate
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = ceil(H/Ts); % Horizon steps
            
            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.3)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar(1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % Define Q and R matrices 
            Q = eye(nx);
            Q(2,2) = 170;
            
            R = eye(nu);
            
            % Define U constraints for (both with vectors and scalars)
            uss = 56.6667;
            P_avg_lim_low = 50 - uss;
            P_avg_lim_high = 80 - uss;
            M = [1;-1];
            m = [P_avg_lim_high;-P_avg_lim_low];
            
            % Compute LQR invariant set and final cost
            [~,P,~] = dlqr(mpc.A,mpc.B,Q,R);
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            obj = 0;
            con = [];
            
            for i=1:N-1
                % Discrete Time model constraint, taking into account the disturbance
                con = [con, X(:,i+1) == mpc.A*X(:,i) + mpc.B*U(:,i)+ mpc.B*d_est];
                % Constraints on U
                con = [con, M*U(:,i) <= m];
                % Increment the objective function
                obj = obj + (X(:,i) - x_ref)'*Q*(X(:,i) - x_ref) + (U(:,i) - u_ref)'*R*(U(:,i) - u_ref); 
            end
            % Increment the objective function with the final cost
            obj = obj + (X(:,N) - x_ref)'*P*(X(:,N) - x_ref);
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref, d_est}, U(:,1));
        end
        
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);
            nu = size(mpc.B,2); 
            
            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar(nu,1);
            
            % Reference position (Ignore this before Todo 3.3)
            ref = sdpvar;
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            uss = 56.6667;
            P_avg_lim_low = 50 - uss;
            P_avg_lim_high = 80 - uss;
            M = [1;-1];
            m = [P_avg_lim_high;-P_avg_lim_low];
            
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            % Implement these constraints
            obj = us'*us;
            % Steady state constraints, taking into account the disturbance
            con = [mpc.A*xs + mpc.B*us + mpc.B*d_est == xs];
            % Reference constraints
            con = [con, mpc.C*xs + mpc.D*us == ref];
            % Input constraints
            con = [con, M*us <= m];
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
        end
        
        
        % Compute augmented system and estimator gain for input disturbance rejection
        function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
            
            %%% Design the matrices A_bar, B_bar, L, and C_bar
            %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
            %%% converges to the correct state and constant input disturbance
            %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            [nx, nu] = size(mpc.B);
            ny   = size(mpc.C,1);
            
            %Building the A_bar, B_bar and C_bar matrices
            A_bar = [mpc.A     , mpc.B;
                    zeros(1,nx),1           ];

            B_bar = [mpc.B;zeros(1,nu)];
            C_bar = [mpc.C,zeros(ny,1)];
            
            %Placing of the poles for the closed loop controller: ( A_bar + L*C_bar )
            L = -place(A_bar',C_bar',[0.1,0.2,0.3])';

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
    end
end
