hold on, grid on

Est_Wind_Height = [W0_vec(start_filter:end,1:2), heights(start_filter:end)];
Est_Wind_Height = sortrows(Est_Wind_Height, size(Est_Wind_Height,2));
x_est = Est_Wind_Height(:,end);
dir_est = atan2d(Est_Wind_Height(:,2),Est_Wind_Height(:,1));
plot(dir_est,x_est, 'o');

% Real(Flown) Wind Points
Act_Wind_Height = [W_log.signals.values(N_start+start_filter:N_end,1:2), heights(N_start+start_filter:end)];
Act_Wind_Height = sortrows(Act_Wind_Height, 3);
x_act = Act_Wind_Height(:,end);
dir_act = atan2d(Act_Wind_Height(:,2),Act_Wind_Height(:,1));
plot(dir_act,x_act, 'o');

% Statistical Points
ref_heights = linspace(min(heights), max(heights), segments+1);
Est_WindDir_meanstd = zeros(size(ref_heights,2)-1,3);
for i = 1:size(ref_heights,2)-1
    h_ind = find(Est_Wind_Height(:,3)>= ref_heights(i) & Est_Wind_Height(:,3)< ref_heights(i+1));
    Est_WindDir_meanstd(i,:) = [mean(dir_est(h_ind,1)), std(dir_est(h_ind,1)), ref_heights(i)];
end
errorbar(Est_WindDir_meanstd(:,1),Est_WindDir_meanstd(:,3), Est_WindDir_meanstd(:,2), 'o-', "horizontal", 'Color', 'black', 'MarkerFaceColor','black', 'LineWidth', 1)

xlabel("Wind Direction [deg]"), ylabel("Height [m]");
%xlim([-15 15]);
ylim([10 30]);
title('Wind Direction Profile');
legend('Estimated', 'Real','Statistical', 'Location','northwest');
