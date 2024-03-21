close all

options = optimoptions('fsolve');
options.MaxFunctionEvaluations = 6000;
options.MaxIterations = 2000;
options.Display = 'off';
options.Algorithm = 'levenberg-marquardt';

N_opt = 2500;

W0 = W_log.signals.values(1,:)';
W_vec = [zeros(N_opt,3)];
x_meas_vec = [zeros(N_opt,6)];
Psi_meas_vec = [zeros(N_opt,1)];
Fline_meas_vec = [zeros(N_opt,3)];

% Noise levels, 0.01 is 1% gaussian noise
%x_noiseLvl = [0.01;0.1;0.1;0.1;0.001;0.1];
x_noiseLvl = 0;
Psi_noiseLvl = 0;
Fline_noiseLvl = 0;

% Generating Parameters

W_meas               = W_log.signals.values(1,:)';
    
x_meas               = states.signals.values(1,:)';
x_meas(1)            = pi/2 - x_meas(1);
x_meas(3)            = -x_meas(3);
x_meas               = x_meas + x_noiseLvl.*x_meas.*randn(size(x_meas)); %adding noise to measurement
Psi_meas             = Psi.signals.values(1);
Psi_meas             = Psi_meas + Psi_noiseLvl*Psi_meas.*randn(size(Psi_meas));
Fline_meas           = [0; 0; Forces.signals.values(1,end)];
Fline_meas           = Fline_meas + Fline_noiseLvl*Fline_meas.*randn(size(Fline_meas));
xdot_meas            = statesdot.signals.values(1,:)';

pm.A                 = area;
pm.m                 = mass;
pm.rho               = rho;
pm.g                 = g;
pm.n_l               = n_line;
pm.d_l               = Line_diameter;
pm.CD_l              = CD_Line;
pm.rho_l             = Line_density;
pm.alpha_0           = alpha_0;
pm.alpha_max         = alpha_max;
pm.alpha_min         = alpha_min;
pm.alpha_var         = alpha_var;
pm.Cl_var            = CL_var;
pm.Cd_var            = CD_var;

codegen kiteEq -args {x_meas,Psi_meas,Fline_meas,W_meas,pm} -lang:c++

tic
for ind = 1:N_opt

    W_meas               = W_log.signals.values(ind,:)';
    
    x_meas               = states.signals.values(ind,:)';
    x_meas(1)            = pi/2 - x_meas(1);
    x_meas(3)            = -x_meas(3);
    x_meas               = x_meas + x_noiseLvl.*x_meas.*randn(size(x_meas)); %adding noise to measurement
    Psi_meas             = Psi.signals.values(ind);
    Psi_meas             = Psi_meas + Psi_noiseLvl*Psi_meas.*randn(size(Psi_meas));
    Fline_meas           = [0; 0; Forces.signals.values(ind,end)];
    Fline_meas           = Fline_meas + Fline_noiseLvl*Fline_meas.*randn(size(Fline_meas));
    xdot_meas            = statesdot.signals.values(ind,:)';

    x_meas_vec(ind,:)    = x_meas';
    Psi_meas_vec(ind,:)  = Psi_meas;
    Fline_meas_vec(ind,:)= Fline_meas';
    
    fun = @(W)accDiff(x_meas,Psi_meas,Fline_meas,W,pm, xdot_meas);
    [W,~,exit,out] = fsolve(fun, W0, options); 
    %exit = 1;
    W_vec(ind,:) = W'; 

    W0 = W;   % soft start

    % Prints
    fprintf("Step %d, Estimated Wind is [%-2.2f,%-2.2f,%-2.2f], Actual Wind is [%-2.2f,%-2.2f,%-2.2f]", ...
        ind-1, W(1), W(2), W(3), W_meas(1),W_meas(2),W_meas(3));

    if exit == 0
        fprintf(" Exceeded, iter: %d, fncCount: %d", out.iterations, out.funcCount);
    end
    if exit == -2 || exit == -3
        fprintf(" Eq not solved, exit: %d, msg: \n %s \n", exit, out.message);
    end
    fprintf("\n");
    
end
toc

RMSE = rmse(W_log.signals.values(1:N_opt,:), W_vec(:,:));
fprintf("RMSE: %f, %f, %f, sum: %f\n", RMSE, sum(RMSE));

Xtime = (1:N_opt)/100;

figure(1)
plot(Xtime, W_log.signals.values(1:N_opt,1), Xtime, W_vec(:,1), 'LineWidth', 1);
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_x$', 'Interpreter','latex');
legend('Actual $W_x$','Estimated $W_x$', 'Interpreter', 'latex');
ylim([6 16]);

figure(2)
plot(Xtime, W_log.signals.values(1:N_opt,2), Xtime, W_vec(:,2), 'LineWidth', 1);
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_y$', 'Interpreter','latex');
legend('Actual $W_y$','Estimated $W_y$', 'Interpreter', 'latex');

% figure(3)
% plot(Xtime, W_log.signals.values(1:N_opt,3),Xtime, W_vec(:,3), 'LineWidth', 1);
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Speed (m/s)', 'Interpreter','latex');
% title('$W_z$', 'Interpreter','latex');
% legend('Actual $W_z$','Estimated $W_z$', 'Interpreter', 'latex', 'Location','southeast');
% h = get(gca,'Children');
% set(gca,'Children',[h(2) h(1)]);
% hold off
% 
% figure(4)
% subplot(3,2,1)
% plot(Xtime, x_meas_vec(1:N_opt,1), Xtime, pi/2 - states.signals.values(1:N_opt,1));
% title('theta')
% subplot(3,2,2)
% plot(Xtime, x_meas_vec(1:N_opt,2), Xtime, states.signals.values(1:N_opt,2));
% title('phi')
% subplot(3,2,3)
% plot(Xtime, x_meas_vec(1:N_opt,3), Xtime, -states.signals.values(1:N_opt,3));
% title('thetadot')
% subplot(3,2,4)
% plot(Xtime, x_meas_vec(1:N_opt,4), Xtime, states.signals.values(1:N_opt,4));
% title('phidot')
% subplot(3,2,5)
% plot(Xtime, x_meas_vec(1:N_opt,5), Xtime, states.signals.values(1:N_opt,5));
% title('L')
% subplot(3,2,6)
% plot(Xtime, x_meas_vec(1:N_opt,6), Xtime, states.signals.values(1:N_opt,6));
% title('Ldot')
% 
% figure(5)
% subplot(2,1,1)
% plot(Xtime, Psi_meas_vec(:), Xtime, Psi.signals.values(1:N_opt));
% title('Psi')
% subplot(2,1,2)
% plot(Xtime, Fline_meas_vec(:,3), Xtime, Forces.signals.values(1:N_opt,end));
% title('Fline')


function F = accDiff(x_meas,Psi_meas,Fline_meas,W,pm, xdot_meas)
    
    xdot = kiteEq_mex(x_meas,Psi_meas,Fline_meas,W,pm);
    F = xdot_meas - xdot;

end
