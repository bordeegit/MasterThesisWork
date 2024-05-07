close all 
%clear
% Used for Hyperparameters optimization, remove clear before
%clearvars -except parameters Cd_mod Cd_mean RMSE_cell W0_cell iter

load FlightData/Standard.mat

% Structural Parameters
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;
parameters.Cd           = mean(Cd.signals.values);
parameters.Cl           = mean(Cl.signals.values);


% Size Initialization for codegen
parameters.r_meas       = Position.signals.values(1,:)';
parameters.rd_meas      = PositionDot.signals.values(1,:)';
parameters.rdd_meas     = PositionDotDot.signals.values(1,:)';
parameters.F_T_norm     = Forces.signals.values(1,end);


% Simulation/Optimization Parameters

N_start                 = 1;
N_opt                   = 5000; % Number of steps to perform optimization
printFlag               = true;
codegenFlag             = true;

z0                      = [14;5;0;%W_log.signals.values(N_start,:)';
                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
%z0                      = [8;3;0; %W_log.signals.values(1,:)';
%                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
Nop                     = size(z0,1);

N_end                   = N_start+N_opt-1;
W0_vec                  = zeros(N_opt,3);
zl_vec                  = zeros(N_opt,3);
heights                 = zeros(N_opt,1);
parameters.Q            = 5e2*diag(ones(3,1));
parameters.Qw           = 1e7*diag(ones(6,1));
%parameters.gamma        = 1e4;
parameters.zold         = z0;

% Linear Inequality Constraints
A = [];
b = [];

% Linear Equality Constaints (W_z = 0)
%Useless to use this kind of constraint, moved it to bound (more efficient)
Aeq = [];
beq = [];

% Bounds on W_x, W_y
lb = [5;-2;0;-Inf;-Inf;-Inf];
ub = [15;10;0;Inf;Inf;Inf];


% Optimization Options
options                             = optimoptions('fmincon');
options.Algorithm                   = 'sqp';
options.FiniteDifferenceType        = 'forward';
options.Display                     = 'off';

options.StepTolerance               = 1e-10;
options.OptimalityTolerance         = 1e-10;
options.ConstraintTolerance         = 1e-10;
options.MaxIterations               = 500;
options.MaxFunctionEvaluations      = 5000; 
options.SpecifyObjectiveGradient    = true;
options.SpecifyConstraintGradient   = true;
options.CheckGradients              = false;


if codegenFlag
    codegen Wind_cost -args {z0,parameters} -lang:c++
    codegen normconstr -args {z0, parameters.rd_meas} -lang:c++
end

% Filtering AoA
%filteredAlpha = medfilt1(alpha.signals.values,10);

tic
for i = N_start:N_end
    parameters.r_meas       = Position.signals.values(i,:)';
    parameters.rd_meas      = PositionDot.signals.values(i,:)';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
    parameters.F_T_norm     = Forces.signals.values(i,end);
    %parameters.Cl           = interp1(alpha_var,CL_var,filteredAlpha(i),'spline','extrap');
    %parameters.Cd           = interp1(alpha_var,CD_var,filteredAlpha(i),'spline','extrap');

    nl_con = @(z)normconstr_mex(z, parameters.rd_meas);
    fun = @(z)Wind_cost_mex(z,parameters);
    [zstar,~,exitflag,out] = fmincon(fun,z0,A,b,Aeq,beq,lb,ub,nl_con,options);
    z0                      = zstar;
    W0_vec(i-N_start+1,:)   = zstar(1:3)';
    heights(i-N_start+1,:)  = parameters.r_meas(3);
    zl_vec(i-N_start+1,:)   = zstar(4:6)';
    parameters.zold         = zstar;
    if printFlag
        fprintf("Iteration %d done, Wind is [%7.4f %7.4f %7.4f], (norm %f), iter: %3d, feval: %3d, exit:%2d \n", ...
            i, zstar(1), zstar(2), zstar(3), norm(zstar(4:6)), out.iterations, out.funcCount, exitflag);
    end
    %fprintf("Cd : %7.4f   Cl : %7.4f\n", parameters.Cl, parameters.Cd);
            
end
toc

RMSE = rmse(W_log.signals.values(N_start:N_end,:), W0_vec(:,:));
fprintf("RMSE: %f, %f, %f, 2-norm: %f\n", RMSE, norm(RMSE));

%% Printing results

Xtime = (N_start:N_end)/100;

% Difference in Wind Magnitude in X direction
figure(1)
plot(Xtime, W_log.signals.values(N_start:N_end,1), '--' , Xtime, W0_vec(:,1));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_x$', 'Interpreter','latex');
legend('Actual $W_x$','Estimated $W_x$', 'Interpreter', 'latex');
ylim([6 16]);
fig1 = gca;

% Difference in Wind Magnitude in Y direction
figure(2)
plot(Xtime, W_log.signals.values(N_start:N_end,2),'--' ,Xtime, W0_vec(:,2));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_y$', 'Interpreter','latex');
legend('Actual $W_y$','Estimated $W_y$', 'Interpreter', 'latex');
fig2 = gca;

% Difference in Norm between estimated and actual wind
figure(3)
plot(Xtime, vecnorm(W_log.signals.values(N_start:N_end,:)'), '--'), hold on
plot(Xtime, vecnorm(W0_vec')), hold off
xlabel('Time (s)','Interpreter','latex');
title('Wind Norm', 'Interpreter','latex');
legend('Actual $|W|$','Estimated $|W|$', 'Interpreter', 'latex');

% 3D Position of trajectory
figure(4),
plot3(Position.signals.values(:,1),Position.signals.values(:,2),Position.signals.values(:,3),'k'),grid on,hold on
plot3(0,0,0,'k*')
plot3([0, Position.signals.values(end,1)],...
     [0, Position.signals.values(end,2)],...
     [0, Position.signals.values(end,3)], 'm-o')
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)'),
%indices = find(W0_vec(:,2) > 3.2 | W0_vec(:,2) < -1.3);
%indices = indices(indices > 2000);
%test = Position.signals.values(indices,:);
%plot3(test(:,1),test(:,2),test(:,3),'or'); hold off
ind_pos = find(W0_vec(:,2) > 3.2);
ind_neg = find(W0_vec(:,2) < -1.3);
ind_pos = ind_pos(ind_pos > 2000);
ind_neg = ind_neg(ind_neg > 2000);
test_pos = Position.signals.values(ind_pos,:);
test_neg = Position.signals.values(ind_neg,:);
plot3(test_pos(:,1),test_pos(:,2),test_pos(:,3),'or'); 
plot3(test_neg(:,1),test_neg(:,2),test_neg(:,3),'ob');
hold off

% Euclidean distance between Lift estimated and actual direction
figure(5);
Fl = Forces.signals.values(N_start:N_opt, 4:6);
zl_distance = zeros(N_opt);
for i = 1:N_opt
    zl_distance(i) = norm(Fl(i,:)/norm(Fl(i,:)) - zl_vec(i));
end
plot(Xtime, zl_distance)
xlabel('Time (s)','Interpreter','latex');
ylabel('Magnitude', 'Interpreter','latex');
title('$z_l\: Euclidean\: Norm $', 'Interpreter','latex');


printWindX;

% zl_vec = Forces.signals.values(N_start:N_end,4:6)/norm(Forces.signals.values(N_start:N_end,4:6));
% figure(5);
% subplot(3,1,1);
% plot(Xtime, Forces.signals.values(N_start:N_end,4), '--' , Xtime, zl_vec(:,1));
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Speed (m/s)', 'Interpreter','latex');
% title('$W_x$', 'Interpreter','latex');
% legend('Actual $W_x$','Estimated $W_x$', 'Interpreter', 'latex');




% figure(3)
% plot(Xtime, W_log.signals.values(N_start:N_end,3),'--' , Xtime, W0_vec(:,3));
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Speed (m/s)', 'Interpreter','latex');
% title('$W_z$', 'Interpreter','latex');
% legend('Actual $W_z$','Estimated $W_z$', 'Interpreter', 'latex', 'Location','southeast');
% h = get(gca,'Children');
% set(gca,'Children',[h(2) h(1)]);
% fig3 = gca;

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
