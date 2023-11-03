% Definition of Hyperparameters
Q_values = [1e3, 1e4, 1e5];
Qw_values = [1e5, 1e6, 1e7];
gamma_values = [1e4, 1e5, 1e6];

best_RMSE = inf;
best_hyperparamRMSE = [];
best_time = inf;
best_hyperparamTime = [];

%Constants
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;
parameters.Cl           = mean(Cl.signals.values);
parameters.Cd           = mean(Cd.signals.values);
parameters.deltaNorm    = 0.001;  % Delta for Zl nonlin. ineq. constrain

N_start                 = 3001;
N_opt                   = 150; % Number of steps to perform optimization
N_end                   = N_start+N_opt-1;

codegen Wind_cost -args {z0,parameters} -lang:c++
codegen normconstr -args {z0,parameters.deltaNorm} -lang:c++

for Q = Q_values
    for Qw = Qw_values
        for gamma = gamma_values
       
            z0                      = [W_log.signals.values(N_start,:)';
                                       0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];

            Nop                     = size(z0,1);
            
            W0_vec                  = [zeros(N_opt,3)];
            parameters.Q            = Q*diag(ones(3,1));
            parameters.Qw           = Qw*diag(ones(6,1));
            parameters.gamma        = gamma;
            parameters.zold         = z0;

            
            % Linear Inequalities
            A = [];
            b = [];
            
            % Equality Constaints (W_z = 0)
            Aeq = zeros(Nop,Nop);
            Aeq(3,3) = 1;
            beq = zeros(Nop,1);
            
            % Bounds on W_x, W_y
            lb = [5;-2;-Inf;-Inf;-Inf;-Inf];
            ub = [15;10;Inf;Inf;Inf;Inf];
            
            
            % Optimization Options
            options               = optimset('fmincon');
            options.FinDiffType   = 'forward';
            %options.FinDiffRelStep= 2^-17;
            options.Display       = 'off';
            options.Algorithm     = 'sqp';
            options.MaxIter       = 100;
            options.MaxFunEvals   = 1000; 
            options.GradConstr    = 'on';
            
            %codegen Wind_cost -args {z0,parameters} -lang:c++
            %codegen normconstr -args {z0,parameters.deltaNorm} -lang:c++
            
            % Nonlinear Constrains
            nonlcon = @(z)normconstr_mex(z,parameters.deltaNorm);
            
            tic
            for i = N_start:N_end
                parameters.r_meas       = Position.signals.values(i,:)';
                parameters.rd_meas      = PositionDot.signals.values(i,:)';
                parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
                parameters.F_T_norm     = Forces.signals.values(i,end);
            
                fun = @(z)Wind_cost_mex(z,parameters);
                [zstar,~,exitflag,out] = fmincon(fun,z0,[],[],Aeq,beq,lb,ub,nonlcon,options);
                z0                      = zstar;
                W0_vec(i-N_start+1,:)             = zstar(1:3)';
                parameters.zold         = zstar;
            end
            t = toc;
            
            
            RMSE = rmse(W_log.signals.values(N_start:N_end,:), W0_vec(:,:));

            fprintf("Iteration Ended, RMSE: %f\n", sum(RMSE));
            if sum(RMSE) < best_RMSE
                best_RMSE = sum(RMSE);
                best_hyperparamRMSE = [Q, Qw, gamma];
            end
            
            if t < best_time
                best_time = t;
                best_hyperparamTime = [Q, Qw, gamma];
            end

        end
    end
end
fprintf("RMSE   Best Hyperparameters: Q=%e, Qw=%e, gamma=%e, RMSE=%f\n", ...
            best_hyperparamRMSE(1), best_hyperparamRMSE(2), best_hyperparamRMSE(3), best_RMSE);
fprintf("TIME   Best Hyperparameters: Q=%e, Qw=%e, gamma=%e, time=%f\n", ...
            best_hyperparamTime(1), best_hyperparamTime(2), best_hyperparamTime(3), t);


