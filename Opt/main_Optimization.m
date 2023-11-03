close all 

% Structural Parameters
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;
parameters.Cl           = mean(Cl.signals.values);
parameters.Cd           = mean(Cd.signals.values);

parameters.rd_meas      = PositionDot.signals.values(1,:)';


% Simulation/Optimization Parameters

N_start                 = 2001;
N_opt                   = 500; % Number of steps to perform optimization

z0                      = [W_log.signals.values(N_start,:)';
                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
%z0                      = [8;3;0; %W_log.signals.values(1,:)';
%                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
Nop                     = size(z0,1);

N_end                   = N_start+N_opt-1;
W0_vec                  = [zeros(N_opt,3)];
parameters.Q            = 1e4*diag(ones(3,1));
parameters.Qw           = 1e6*diag(ones(6,1));
parameters.gamma        = 1e1;
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
options                         = optimoptions('fmincon');
options.Algorithm               = 'sqp';
options.FiniteDifferenceType    = 'forward';
options.Display                 = 'off';

options.StepTolerance           = 1e-8;
options.ConstraintTolerance     = 1e-8;
options.OptimalityTolerance     = 1e-8;
options.ConstraintTolerance     = 1e-8;
options.MaxIterations           = 500;
options.MaxFunctionEvaluations  = 5000; 
options.GradConstr              = 'on';

codegen Wind_cost -args {z0,parameters} -lang:c++
codegen normconstr -args {z0, parameters.rd_meas} -lang:c++

% Nonlinear Constrains


tic
for i = N_start:N_end
    parameters.r_meas       = Position.signals.values(i,:)';
    parameters.rd_meas      = PositionDot.signals.values(i,:)';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
    parameters.F_T_norm     = Forces.signals.values(i,end);

    nonlcon = @(z)normconstr_mex(z, parameters.rd_meas);
    fun = @(z)Wind_cost_mex(z,parameters);
    [zstar,~,exitflag,out] = fmincon(fun,z0,[],[],Aeq,beq,lb,ub,nonlcon,options);
    z0                      = zstar;
    W0_vec(i-N_start+1,:)   = zstar(1:3)';
    parameters.zold         = zstar;
    fprintf("Iteration %d done, Wind is [%f %f %f], (norm %f), iter: %d, feval: %d, exit:%d \n", ...
            i, zstar(1), zstar(2), zstar(3), norm([zstar(4) zstar(5) zstar(6)]), out.iterations, out.funcCount, exitflag);
end
toc


RMSE = rmse(W_log.signals.values(N_start:N_end,:), W0_vec(:,:));
fprintf("RMSE: %f, %f, %f, sum: %f\n", RMSE, sum(RMSE));

figure(1)
plot(N_start:N_end, W_log.signals.values(N_start:N_end,1), N_start:N_end, W0_vec(:,1));
figure(2)
plot(N_start:N_end, W_log.signals.values(N_start:N_end,2), N_start:N_end, W0_vec(:,2));
%figure(3)
%plot(N_start:N_end, W_log.signals.values(N_start:N_end,3), N_start:N_end, W0_vec(:,3));
