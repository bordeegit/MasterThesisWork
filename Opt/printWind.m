Xtime = (N_start:N_end)*T_s;
f = figure;
f.Position = [300 300 1200 500];
% Difference in Wind Magnitude in X direction
subplot(1,2,1)
plot(Xtime, W(N_start:N_end,1), '--' , Xtime, W0_vec(:,1));
xlabel('Time (s)','Interpreter','latex'), grid on;
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_x$', 'Interpreter','latex');
legend('Actual $W_x$','Estimated $W_x$', 'Interpreter', 'latex');
grid on %,ylim([6 16]);

% Difference in Wind Magnitude in Y direction
subplot(1,2,2)
plot(Xtime, W(N_start:N_end,2),'--' ,Xtime, W0_vec(:,2));
xlabel('Time (s)','Interpreter','latex');
ylabel('Speed (m/s)', 'Interpreter','latex');
title('$W_y$', 'Interpreter','latex');
legend('Actual $W_y$','Estimated $W_y$', 'Interpreter', 'latex');
ylim([-3 11]), grid on;