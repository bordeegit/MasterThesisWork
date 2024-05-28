close all 

set(groot,'DefaultAxesFontSize', 15);
set(groot,'DefaultLineLineWidth', 1.5);
set(groot,'DefaultTextInterpreter', 'Latex');
set(groot,'DefaultAxesTickLabelInterpreter', 'Latex');  
set(groot,'DefaultLegendInterpreter', 'Latex');

%clear
% Used for Hyperparameters optimization, remove clear before
%clearvars -except parameters Cd_mod Cd_mean RMSE_cell W0_cell iter

%load FlightData/Standard_LinY.mat

% Structural Parameters
parameters.rho          = rho;
parameters.A            = area;
parameters.mk           = mass;
parameters.mt_noL       = mt_noL;
parameters.g            = g;
parameters.Cd_l         = CD_Line;
parameters.d_l          = Line_diameter;
parameters.n_l          = n_line;
parameters.Cd           = mean(Cd.signals.values);
parameters.Cl           = mean(Cl.signals.values);


% Size Initialization for codegen
parameters.r_meas       = Position.signals.values(1,:)';
parameters.rd_meas      = PositionDot.signals.values(1,:)';
parameters.rdd_meas     = PositionDotDot.signals.values(1,:)';
parameters.F_T_norm     = Forces.signals.values(1,end);


% Simulation/Optimization Parameters

N_start                 = 1;
N_opt                   = 2500; % Number of steps to perform optimization
printFlag               = true;
codegenFlag             = false;

z0                      = [6;0;%W_log.signals.values(N_start,1:2)';
                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)];

Nop                     = size(z0,1);
N_end                   = N_start+N_opt-1;
W0_vec                  = zeros(N_opt,3);
zl_vec                  = zeros(N_opt,3);
heights                 = zeros(N_opt,1);
parameters.Q            = 5e2*diag(ones(3,1)); %5e2
parameters.Qw           = 1e7*diag(ones(5,1)); %1e7
parameters.zold         = z0;

% Linear Equality/Inequality Constraints
Aeq = [];
beq = [];
A = [];
b = [];

% Bounds on W_x, W_y (can also add bound on components of zl)
lb = [5;-2;-1;-1;-1];
ub = [15;20;1;1;1];


% Optimization Options
options                             = optimoptions('fmincon');
options.Algorithm                   = 'sqp';
options.FiniteDifferenceType        = 'central';
options.Display                     = 'off';

% options.StepTolerance               = 1e-10;
% options.OptimalityTolerance         = 1e-10;
% options.ConstraintTolerance         = 1e-10;
% options.MaxIterations               = 500;
% options.MaxFunctionEvaluations      = 5000; 
options.SpecifyObjectiveGradient    = true;
options.SpecifyConstraintGradient   = true;


if codegenFlag
    codegen Wind_cost -args {z0,parameters} -lang:c++
    codegen normconstr -args {z0, parameters.rd_meas} -lang:c++
end

% Check Derivatives (v2024 only)
if isMATLABReleaseOlderThan("R2024a")
    if options.SpecifyObjectiveGradient
        assert(checkGradients(@(z)Wind_cost(z, parameters), randn(5,1), options, Display="on"))
    end
    if options.SpecifyConstraintGradient
        assert(all(checkGradients(@(z)normconstr(z, parameters.rd_meas), randn(5,1), options, IsConstraint=true, Display="on")))
    end
end

%Add noise to Measurements (0.01 = 1%)
noiseLvl = 0.1;

tic
for i = N_start:N_end
    parameters.r_meas       = Position.signals.values(i,:)' + noiseLvl.*Position.signals.values(i,:)'.*randn(size(Position.signals.values(i,:)))';
    parameters.rd_meas      = PositionDot.signals.values(i,:)' + noiseLvl.*PositionDot.signals.values(i,:)'.*randn(size(PositionDot.signals.values(i,:)))';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)' + noiseLvl.*PositionDotDot.signals.values(i,:)'.*randn(size(PositionDotDot.signals.values(i,:)))';
    parameters.F_T_norm     = Forces.signals.values(i,end) + noiseLvl*Forces.signals.values(i,end)*randn(size(Forces.signals.values(i,end)));

    nl_con = @(z)normconstr_mex(z, parameters.rd_meas);
    fun = @(z)Wind_cost_mex(z,parameters);
    [zstar,~,exitflag,out] = fmincon(fun,z0,A,b,Aeq,beq,lb,ub,nl_con,options);
    z0                      = zstar;
    W0_vec(i-N_start+1,:)   = [zstar(1:2)' 0];
    heights(i-N_start+1,:)  = Position.signals.values(i,3);
    zl_vec(i-N_start+1,:)   = zstar(3:5)';
    parameters.zold         = zstar;
    if printFlag
        fprintf("Iteration %4d done, EstWind is [%7.4f %7.4f], (norm %4.3f), iter: %3d, feval: %3d, exit:%2d \n", ...
            i, zstar(1), zstar(2), norm(zstar(3:5)), out.iterations, out.funcCount, exitflag);
    end     
end
toc

RMSE = rmse(W_log.signals.values(N_start:N_end,:), W0_vec(:,:));
fprintf("RMSE: %f, %f, %f, 2-norm: %f\n", RMSE, norm(RMSE));

%% Printing results

%%% Wind Estimation over time
printWind

%%% 3D Position of trajectory
printTraj

%%% Wind Profile 
printWindProfile

%%% Difference in Norm between estimated and actual wind
% figure(3)
% plot(Xtime, vecnorm(W_log.signals.values(N_start:N_end,:)'), '--'), hold on
% plot(Xtime, vecnorm(W0_vec')), hold off
% xlabel('Time (s)','Interpreter','latex');
% title('Wind Norm', 'Interpreter','latex');
% legend('Actual $|W|$','Estimated $|W|$', 'Interpreter', 'latex');

%%% Finding where peaks are along the trajectory 
% ind_pos = find(W0_vec(:,2) > 3.2);
% ind_neg = find(W0_vec(:,2) < -1.3);   
% ind_pos = ind_pos(ind_pos > 2000);
% ind_neg = ind_neg(ind_neg > 2000);
% test_pos = Position.signals.values(ind_pos,:);
% test_neg = Position.signals.values(ind_neg,:);
% plot3(test_pos(:,1),test_pos(:,2),test_pos(:,3),'or'); 
% plot3(test_neg(:,1),test_neg(:,2),test_neg(:,3),'ob');
% hold off

%%% Euclidean distance between Lift estimated and actual direction
% figure;
% Fl = Forces.signals.values(N_start:N_end, 4:6);
% zl_dot = zeros(N_opt,1);
% for i = N_start:N_end
%     zl_dot(i-N_start+1) = dot(Fl(i-N_start+1,:)/norm(Fl(i-N_start+1,:)),zl_vec(i-N_start+1,:));
% end
% plot(Xtime, zl_dot)
% xlabel('Time (s)','Interpreter','latex'); 
% ylabel('Magnitude', 'Interpreter','latex'), ylim([0.95 1.05]);
% title('$z_l\: Dot\: Product $', 'Interpreter','latex');

%%% Aerodynamic coefficients Hyperparamter Optimization
% Xtime = (N_start:N_end)/100;
% 
% meanCd = parameters.Cd;
% meanCl = parameters.Cl;
%
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
%
% exportgraphics(gca,'wx.pdf','ContentType','vector')
