clc
clear
close all

Cl_mod = [1 0.95 0.90 0.85 0.80 1.05 1.10 1.15 1.20];
Cl_mean = 0.739382164633603; 

RMSE_cell = cell(1, size(Cl_mod,1));
W0_cell = cell(1, size(Cl_mod,1));

iter = 0; 

for Cl_mul = Cl_mod
    iter = iter + 1;
    parameters.Cl = Cl_mul*Cl_mean;
    main_Optimization
    RMSE_cell{iter} = RMSE;
    W0_cell{iter} = W0_vec;
end

%% Plotting
Xtime = (N_start:N_end)/100;

%% CL
%% Wx
figure(1)
plot(Xtime, W_log.signals.values(N_start:N_end,1), 'LineWidth', 1.5), hold on
plot(Xtime, W0_cell{1}(:,1));
for i = 2:5
    plot(Xtime, W0_cell{i}(:,1));
end
title('Longitudinal Wind, Decreasing $C_L$', 'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
legend('Measured $W_x$','$\overline{C_L}$','-5\% $\overline{C_L}$','-10\% $\overline{C_L}$',...
    '-15\% $\overline{C_L}$','-20\% $\overline{C_L}$', 'Interpreter', 'latex', 'Location', 'none', ...
    'Position', [0.168332885077784 0.716067207757645 0.106734229844432 0.190998917818044]);
set(gcf,'Position',[100 100 1000 500])
exportgraphics(gca,'WxCl_down.pdf','ContentType','vector')
hold off 


figure(2)
plot(Xtime, W_log.signals.values(N_start:N_end,1), 'LineWidth', 1.5), hold on
plot(Xtime, W0_cell{1}(:,1));
for i = 6:9
    plot(Xtime, W0_cell{i}(:,1));
end
title('Longitudinal Wind, Increasing $C_L$', 'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
legend('Measured $W_x$','$\overline{C_L}$','+5\% $\overline{C_L}$','+10\% $\overline{C_L}$',...
    '+15\% $\overline{C_L}$','+20\% $\overline{C_L}$', 'Interpreter', 'latex', 'Location', 'none', ...
    'Position', [0.168332885077784 0.716067207757645 0.106734229844432 0.190998917818044]);
set(gcf,'Position',[100 100 1000 500])
exportgraphics(gca,'WxCl_up.pdf','ContentType','vector')
hold off 

%%% Wy
figure(3)
plot(Xtime, W_log.signals.values(N_start:N_end,2), 'LineWidth', 1.5), hold on
plot(Xtime, W0_cell{1}(:,2));
for i = 2:5
    plot(Xtime, W0_cell{i}(:,2));
end
title('Lateral Wind, Decreasing $C_L$', 'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
legend('Measured $W_y$','$\overline{C_L}$','-5\% $\overline{C_L}$','-10\% $\overline{C_L}$',...
    '-15\% $\overline{C_L}$','-20\% $\overline{C_L}$', 'Interpreter', 'latex', 'Location', 'northeast');
set(gcf,'Position',[100 100 1000 500])
exportgraphics(gca,'WyCl_down.pdf','ContentType','vector')
hold off 

figure(4)
plot(Xtime, W_log.signals.values(N_start:N_end,2), 'LineWidth', 1.5), hold on
plot(Xtime, W0_cell{1}(:,2));
for i = 6:9
    plot(Xtime, W0_cell{i}(:,2));
end
title('Lateral Wind, Increasing $C_L$', 'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
legend('Measured $W_y$','$\overline{C_L}$','+5\% $\overline{C_L}$','+10\% $\overline{C_L}$',...
    '+15\% $\overline{C_L}$','+20\% $\overline{C_L}$', 'Interpreter', 'latex', 'Location', 'northeast');
set(gcf,'Position',[100 100 1000 500])
exportgraphics(gca,'WyCl_up.pdf','ContentType','vector')
hold off 


%RMSE_cell{:}
%cellfun(@(x) sum(x), RMSE_cell)'

% RMSE_mat = zeros(iter,3);
% mapeMate = zeros(iter,3);
% for i = 1:iter
%     RMSE_mat(i,:) = RMSE_cell{i};
%     mapeMat(i,:) = mape(W0_cell{i}, W_log.signals.values(N_start:N_end,:));
% end
% 
% mapeMat
% 
% cellfun(@(x) norm(x), RMSE_cell)'
% 
% mean(RMSE_mat)
% std(RMSE_mat)

figure(1)
plot(Xtime, W_log.signals.values(N_start:N_end,1), 'LineWidth', 1.5), hold on
plot(Xtime, W0_cell{1}(:,1));
plot(Xtime, W0_cell{3}(:,1));
plot(Xtime, W0_cell{7}(:,1));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('Longitudinal Wind, Varying $C_L$', 'Interpreter','latex');
legend('Measured $W_x$','$\overline{C_L}$','-10\% $\overline{C_L}$',...
    '+10\% $\overline{C_L}$', 'Interpreter', 'latex', 'Location', 'north');


