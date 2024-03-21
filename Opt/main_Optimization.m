close all 
clearvars -except parameters Cd_mod Cd_mean RMSE_cell W0_cell iter

load FlightData.mat

% Structural Parameters
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;
parameters.Cd           = mean(Cd.signals.values);
parameters.Cl           = mean(Cl.signals.values);


% Need these for codegen
parameters.r_meas       = Position.signals.values(1,:)';
parameters.rd_meas      = PositionDot.signals.values(1,:)';
parameters.rdd_meas     = PositionDotDot.signals.values(1,:)';
parameters.F_T_norm     = Forces.signals.values(1,end);


% Simulation/Optimization Parameters

N_start                 = 1;
N_opt                   = 2500; % Number of steps to perform optimization
printFlag               = true;
codegenFlag             = false;

z0                      = [W_log.signals.values(N_start,:)';
                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
%z0                      = [8;3;0; %W_log.signals.values(1,:)';
%                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
Nop                     = size(z0,1);

N_end                   = N_start+N_opt-1;
W0_vec                  = [zeros(N_opt,3)];
parameters.Q            = 5e2*diag(ones(3,1));
parameters.Qw           = 1e7*diag(ones(6,1));
%parameters.gamma        = 1e4;
parameters.zold         = z0;

% Linear Inequality Constraints
A = [];
b = [];

% Linear Equality Constaints (W_z = 0)
Aeq = zeros(Nop,Nop);
Aeq(3,3) = 1;
beq = zeros(Nop,1);

% Bounds on W_x, W_y
lb = [5;-2;-Inf;-Inf;-Inf;-Inf];
ub = [15;10;Inf;Inf;Inf;Inf];


% Optimization Options
options                             = optimoptions('fmincon');
options.Algorithm                   = 'sqp';
options.FiniteDifferenceType        = 'forward';
options.Display                     = 'iter';

options.StepTolerance               = 1e-10;
options.ConstraintTolerance         = 1e-10;
options.OptimalityTolerance         = 1e-10;
options.ConstraintTolerance         = 1e-10;
options.MaxIterations               = 500;
options.MaxFunctionEvaluations      = 5000; 
options.SpecifyObjectiveGradient    = true;
options.SpecifyConstraintGradient   = true;
%options.CheckGradients              = true;

if codegenFlag
    codegen Wind_cost -args {z0,parameters} -lang:c++
    codegen normconstr -args {z0, parameters.rd_meas} -lang:c++
end


tic
for i = N_start:N_end
    parameters.r_meas       = Position.signals.values(i,:)';
    parameters.rd_meas      = PositionDot.signals.values(i,:)';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
    parameters.F_T_norm     = Forces.signals.values(i,end);
    %parameters.Cl           = Cl.signals.values(i);
    %parameters.Cd           = Cd.signals.values(i);

    nonlcon = @(z)normconstr_mex(z, parameters.rd_meas);
    fun = @(z)Wind_cost_mex(z,parameters);
    [zstar,~,exitflag,out] = fmincon(fun,z0,A,b,Aeq,beq,[],[],nonlcon,options);
    z0                      = zstar;
    W0_vec(i-N_start+1,:)   = zstar(1:3)';
    parameters.zold         = zstar;
    if printFlag
        fprintf("Iteration %d done, Wind is [%f %f %f], (norm %f), iter: %d, feval: %d, exit:%d \n", ...
            i, zstar(1), zstar(2), zstar(3), norm([zstar(4) zstar(5) zstar(6)]), out.iterations, out.funcCount, exitflag);
    end
            
end
toc


RMSE = rmse(W_log.signals.values(N_start:N_end,:), W0_vec(:,:));
fprintf("RMSE: %f, %f, %f, sum: %f\n", RMSE, sum(RMSE));

Xtime = (N_start:N_end)/100;

figure(1)
plot(Xtime, W_log.signals.values(N_start:N_end,1), Xtime, W0_vec(:,1));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_x$', 'Interpreter','latex');
legend('Actual $W_x$','Estimated $W_x$', 'Interpreter', 'latex');
ylim([6 16]);
fig1 = gca;

figure(2)
plot(Xtime, W_log.signals.values(N_start:N_end,2), Xtime, W0_vec(:,2));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_y$', 'Interpreter','latex');
legend('Actual $W_y$','Estimated $W_y$', 'Interpreter', 'latex');
fig2 = gca;

figure(3)
plot(Xtime, W_log.signals.values(N_start:N_end,3), Xtime, W0_vec(:,3));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_z$', 'Interpreter','latex');
legend('Actual $W_z$','Estimated $W_z$', 'Interpreter', 'latex', 'Location','southeast');
h = get(gca,'Children');
set(gca,'Children',[h(2) h(1)]);
fig3 = gca;

Xtime = (N_start:N_end)/100;

meanCd = parameters.Cd;
meanCl = parameters.Cl;

% figure(1)
% plot(Xtime, Cd.signals.values(N_start:N_end)), hold on
% plot(Xtime, meanCd*ones(size(Xtime)))
% plot(Xtime, 1.20*meanCd*ones(size(Xtime)))
% plot(Xtime, 0.80*meanCd*ones(size(Xtime)))
% legend('$C_D$', '$\overline{C_D}$', '$+20\% \overline{C_D}$', '$-20\% \overline{C_D}$','Interpreter', 'latex');
% xlabel('Time (s)'); 
% exportgraphics(gca,'Cdvar.pdf','ContentType','vector');
% hold off 
% 
% figure(2)
% plot(Xtime, Cl.signals.values(N_start:N_end)), hold on
% plot(Xtime, meanCl*ones(size(Xtime)))
% plot(Xtime, 1.20*meanCl*ones(size(Xtime)))
% plot(Xtime, 0.80*meanCl*ones(size(Xtime)))
% legend('$C_L$', '$\overline{C_L}$', '$+20\% \overline{C_L}$', '$-20\% \overline{C_L}$','Interpreter', 'latex');
% xlabel('Time (s)');
% exportgraphics(gca,'Clvar.pdf','ContentType','vector');
% hold off


%exportgraphics(gca,'wx.pdf','ContentType','vector')
