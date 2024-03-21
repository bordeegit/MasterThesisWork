clc
clear
close all

Cd_mod = [1 0.95 0.90 0.85 0.80 1.05 1.10 1.15 1.20];
Cd_mean = 0.108648114974973;

RMSE_cell = cell(1, size(Cd_mod,1));
W0_cell = cell(1, size(Cd_mod,1));

iter = 0; 

for Cd_mul = Cd_mod
    iter = iter + 1;
    parameters.Cd = Cd_mul*Cd_mean;
    main_Optimization
    RMSE_cell{iter} = RMSE;
    W0_cell{iter} = W0_vec;
end

%% Plotting
Xtime = (N_start:N_end)/100;

%% CD
%%% Wx
% figure(1)
% plot(Xtime, W_log.signals.values(N_start:N_end,1), 'LineWidth', 1.5), hold on
% plot(Xtime, W0_cell{1}(:,1));
% for i = 2:5
%     plot(Xtime, W0_cell{i}(:,1));
% end
% title('Longitudinal Wind, Decreasing $C_D$', 'Interpreter','latex');
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Speed (m/s)', 'Interpreter','latex');
% legend('Measured $W_x$','$\overline{C_D}$','-5\% $\overline{C_D}$','-10\% $\overline{C_D}$',...
%     '-15\% $\overline{C_D}$','-20\% $\overline{C_D}$', 'Interpreter', 'latex', 'Location', 'none', ...
%     'Position', [0.168332885077784 0.716067207757645 0.106734229844432 0.190998917818044]);
% set(gcf,'Position',[100 100 1000 500])
% exportgraphics(gca,'WxCd_down.pdf','ContentType','vector')
% hold off 
% 
% 
% figure(2)
% plot(Xtime, W_log.signals.values(N_start:N_end,1), 'LineWidth', 1.5), hold on
% plot(Xtime, W0_cell{1}(:,1));
% for i = 6:9
%     plot(Xtime, W0_cell{i}(:,1));
% end
% title('Longitudinal Wind, Increasing $C_D$', 'Interpreter','latex');
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Speed (m/s)', 'Interpreter','latex');
% legend('Measured $W_x$','$\overline{C_D}$','+5\% $\overline{C_D}$','+10\% $\overline{C_D}$',...
%     '+15\% $\overline{C_D}$','+20\% $\overline{C_D}$', 'Interpreter', 'latex', 'Location', 'none', ...
%     'Position', [0.168332885077784 0.716067207757645 0.106734229844432 0.190998917818044]);
% set(gcf,'Position',[100 100 1000 500])
% exportgraphics(gca,'WxCd_up.pdf','ContentType','vector')
% hold off 
% 
% %%% Wy
% figure(3)
% plot(Xtime, W_log.signals.values(N_start:N_end,2), 'LineWidth', 1.5), hold on
% plot(Xtime, W0_cell{1}(:,2));
% for i = 2:5
%     plot(Xtime, W0_cell{i}(:,2));
% end
% title('Lateral Wind, Decreasing $C_D$', 'Interpreter','latex');
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Speed (m/s)', 'Interpreter','latex');
% legend('Measured $W_y$','$\overline{C_D}$','-5\% $\overline{C_D}$','-10\% $\overline{C_D}$',...
%     '-15\% $\overline{C_D}$','-20\% $\overline{C_D}$', 'Interpreter', 'latex', 'Location', 'northeast');
% set(gcf,'Position',[100 100 1000 500])
% exportgraphics(gca,'WyCd_down.pdf','ContentType','vector')
% hold off 
% 
% figure(4)
% plot(Xtime, W_log.signals.values(N_start:N_end,2), 'LineWidth', 1.5), hold on
% plot(Xtime, W0_cell{1}(:,2));
% for i = 6:9
%     plot(Xtime, W0_cell{i}(:,2));
% end
% title('Lateral Wind, Increasing $C_D$', 'Interpreter','latex');
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Speed (m/s)', 'Interpreter','latex');
% legend('Measured $W_y$','$\overline{C_D}$','+5\% $\overline{C_D}$','+10\% $\overline{C_D}$',...
%     '+15\% $\overline{C_D}$','+20\% $\overline{C_D}$', 'Interpreter', 'latex', 'Location', 'northeast');
% set(gcf,'Position',[100 100 1000 500])
% exportgraphics(gca,'WyCd_up.pdf','ContentType','vector')
% hold off 


%RMSE_cell{:}
%cellfun(@(x) sum(x), RMSE_cell)

% RMSE_mat = zeros(9,3);
% mapeMate = zeros(iter,3);
% for i = 1:iter
%     RMSE_mat(i,:) = RMSE_cell{i};
%     mapeMat(i,:) = mape(W0_cell{i}, W_log.signals.values(N_start:N_end,:));
% end
% 
% disp(mapeMat)
% cellfun(@(x) norm(x), RMSE_cell)'

figure(1)
plot(Xtime, W_log.signals.values(N_start:N_end,1), 'LineWidth', 1.5), hold on
plot(Xtime, W0_cell{1}(:,1));
plot(Xtime, W0_cell{3}(:,1));
plot(Xtime, W0_cell{7}(:,1));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('Longitudinal Wind, Varying $C_D$', 'Interpreter','latex');
legend('Measured $W_x$','$\overline{C_D}$','-10\% $\overline{C_D}$',...
    '+10\% $\overline{C_D}$', 'Interpreter', 'latex', 'Location', 'north');

