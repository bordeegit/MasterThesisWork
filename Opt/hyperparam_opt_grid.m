% Definition of Hyperparameters
Q_values = [1e1, 5e1, 1e2, 5e2, 5e3];
Qw_values = [1e4, 1e5, 1e6, 1e7, 1e8];
%gamma_values = [1e3, 1e4, 1e5, 1e6, 1e7];

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

N_start                 = 1;
N_opt                   = 3000; % Number of steps to perform optimization
N_end                   = N_start+N_opt-1;

% ONLY for Codegen
z0                      = [W_log.signals.values(N_start,:)';
                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];

Nop                     = size(z0,1);

W0_vec                  = [zeros(N_opt,3)];
parameters.Q            = 1e5*diag(ones(3,1));
parameters.Qw           = 1e5*diag(ones(6,1));
%parameters.gamma        = 1e5;
parameters.zold         = z0;
parameters.r_meas       = Position.signals.values(1,:)';
parameters.rd_meas      = PositionDot.signals.values(1,:)';
parameters.rdd_meas     = PositionDotDot.signals.values(1,:)';
parameters.F_T_norm     = Forces.signals.values(1,end);

codegen Wind_cost -args {z0,parameters} -lang:c++
codegen normconstr -args {z0,parameters.rd_meas} -lang:c++

% Optimization Options
options                         = optimoptions('fmincon');
options.Algorithm               = 'sqp';
options.FiniteDifferenceType    = 'forward';
options.Display                 = 'off';

options.StepTolerance           = 1e-10;
options.ConstraintTolerance     = 1e-10;
options.OptimalityTolerance     = 1e-10;
options.ConstraintTolerance     = 1e-10;
options.MaxIterations           = 500;
options.MaxFunctionEvaluations  = 5000; 
options.SpecifyConstraintGradient   = true;
options.SpecifyObjectiveGradient    = true;

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


for Q = Q_values
    for Qw = Qw_values
       
        z0                      = [W_log.signals.values(N_start,:)';
                                   0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];

        Nop                     = size(z0,1);
        
        W0_vec                  = [zeros(N_opt,3)];
        parameters.Q            = Q*diag(ones(3,1));
        parameters.Qw           = Qw*diag(ones(6,1));
        %parameters.gamma        = gamma;
        parameters.zold         = z0;
        
        tic
        for i = N_start:N_end
            parameters.r_meas       = Position.signals.values(i,:)';
            parameters.rd_meas      = PositionDot.signals.values(i,:)';
            parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
            parameters.F_T_norm     = Forces.signals.values(i,end);
        
            nonlcon = @(z)normconstr_mex(z,parameters.rd_meas);
            fun = @(z)Wind_cost_mex(z,parameters);
            [zstar,~,exitflag,out] = fmincon(fun,z0,A,b,Aeq,beq,lb,ub,nonlcon,options);
            z0                      = zstar;
            W0_vec(i-N_start+1,:)   = zstar(1:3)';
            parameters.zold         = zstar;
        end
        t = toc;
        
        
        RMSE = rmse(W_log.signals.values(N_start:N_end,:), W0_vec(:,:));

        fprintf("Iteration Ended, Q=%e, Qw=%e, RMSE: %f, %f, %f, sum: %f\n",...
                 Q, Qw, RMSE, sum(RMSE));
        if sum(RMSE) < best_RMSE || (sum(RMSE) == best_RMSE && t< best_time)
            best_RMSE = sum(RMSE);
            best_hyperparamRMSE = [Q, Qw];
        end
        
        if t < best_time
            best_time = t;
            best_hyperparamTime = [Q, Qw];
        end

    end
end
fprintf("RMSE   Best Hyperparameters: Q=%e, Qw=%e, RMSE=%f\n", ...
            best_hyperparamRMSE(1), best_hyperparamRMSE(2), best_RMSE);
fprintf("TIME   Best Hyperparameters: Q=%e, Qw=%e, time=%f\n", ...
            best_hyperparamTime(1), best_hyperparamTime(2), t);


