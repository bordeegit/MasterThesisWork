close all 

% Structural Parameters
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;
parameters.Cl           = mean(Cl.signals.values);
parameters.Cd           = mean(Cd.signals.values);


% Simulation/Optimization Parameters

z0                      = [W_log.signals.values(1,:)';
                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
%z0                      = [8;3;0; %W_log.signals.values(1,:)';
%                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
Nop                     = size(z0,1);
N_opt                   = 500; % Number of steps to perform optimization
W0_vec                  = [zeros(N_opt,3)];
parameters.Q            = 1e4*diag(ones(3,1));
parameters.Qw           = 1e6*diag(ones(6,1));
parameters.gamma        = 1e7;
parameters.zold         = z0;
parameters.deltaNorm    = 0.001;  % Delta for Zl nonlin. ineq. constrain

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
options.MaxIter       = 500;
options.MaxFunEvals   = 2000; 

%codegen Wind_cost -args {z0,parameters} -lang:c++
%codegen normconstr -args {z0,parameters.deltaNorm} -lang:c++

% Nonlinear Constrains
nonlcon = @(z)normconstr_mex(z,parameters.deltaNorm);

tic
for i=1:N_opt %i = 2000:2500 %i=1:N_opt
    parameters.r_meas       = Position.signals.values(i,:)';
    parameters.rd_meas      = PositionDot.signals.values(i,:)';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
    parameters.F_T_norm     = Forces.signals.values(i,end);

    fun = @(z)Wind_cost_mex(z,parameters);
    [zstar,~,exitflag,out] = fmincon(fun,z0,[],[],Aeq,beq,lb,ub,nonlcon,options);
    z0                      = zstar;
    W0_vec(i,:)             = zstar(1:3)';
    parameters.zold         = zstar;
    fprintf("Iteration %d done, Wind is [%f %f %f], (norm %f), iter: %d, feval: %d, exit:%d \n", ...
            i, zstar(1), zstar(2), zstar(3), norm([zstar(4) zstar(5) zstar(6)]), out.iterations, out.funcCount, exitflag);
end
toc

RMSE = rmse(W_log.signals.values(1:N_opt,:), W0_vec(:,:));
fprintf("RMSE: %f, %f, %f, sum: %f\n", RMSE, sum(RMSE));

figure(1)
plot(1:N_opt, W_log.signals.values(1:N_opt,1), 1:N_opt, W0_vec(:,1));
figure(2)
plot(1:N_opt, W_log.signals.values(1:N_opt,2), 1:N_opt, W0_vec(:,2));
%figure(3)
%plot(1:N_opt, W_log.signals.values(1:N_opt,3), 1:N_opt, W0_vec(:,3));
