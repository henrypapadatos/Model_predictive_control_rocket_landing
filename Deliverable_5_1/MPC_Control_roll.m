classdef MPC_Control_roll < MPC_Control
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = ceil(H/Ts); % Horizon steps
            
            [nx, nu] = size(mpc.B);
            
            % Steady-state targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % Define Q and R matrices 
            Q = eye(nx);
            Q(1,1) = 100;
            Q(2,2) = 5000;
            
            R = eye(nu);
            
            % Define U constraints for (both with vectors and scalars)
            P_diff_lim = 20;
            M = [1;-1];
            m = [P_diff_lim;P_diff_lim];
            
            % Compute LQR invariant set and final cost
            [~,P,~] = dlqr(mpc.A,mpc.B,Q,R);

            % Question to ask to TA
            % X_lqr = Polyhedron([H; M*K],[h; m]);
            
%             plot_invset(X_f,"Maximum invariant set for the 'sys roll' system");
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            obj = 0;
            con = [];
            
            for i=1:N-1
                % Discrete Time model constraint
                con = [con, X(:,i+1) == mpc.A*X(:,i) + mpc.B*U(:,i)];
                % Constraints on U
                con = [con, M*U(:,i) <= m];
                % Increment the objective function
                % We want to minimize Delta X and Delta U as we do offset
                % tracking
                obj = obj + (X(:,i) - x_ref)'*Q*(X(:,i) - x_ref) + (U(:,i) - u_ref)'*R*(U(:,i) - u_ref); 
            end
%           Increment the objective function with the final cost
            obj = obj + (X(:,N) - x_ref)'*P*(X(:,N) - x_ref);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref}, U(:,1));
        end
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Steady-state targets
            nx = size(mpc.A, 1);
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            % Define U constraints for (both with vectors and scalars)
            P_diff_lim = 20;
            M = [1;-1];
            m = [P_diff_lim;P_diff_lim];
            
            % Implement these constraints
            obj = us'*us;
            % Steady state constraints
            con = [mpc.A*xs + mpc.B*us == xs];
            % Reference constraints
            con = [con, mpc.C*xs + mpc.D*us == ref];
            % Input constraints
            con = [con,M*us <= m];
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end
