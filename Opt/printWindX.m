grid on, hold on

% Input Wind Profile, Interpolations
y_table = 0:0.25:50;
plot(w_x_height_Data*z_scale,height_Data, 'o', 'Color',[0 0 0.5], 'MarkerFaceColor',[0 0 0.5]);
x_table = spline(height_Data,w_x_height_Data*z_scale, y_table);
plot(x_table,y_table, 'Color', [1 0.6 0.2]);

% Estimated Wind Points
Est_Wind_Height = [W0_vec(start_filter:end,1:2), heights(start_filter:end)];
Est_Wind_Height = sortrows(Est_Wind_Height, size(Est_Wind_Height,2));
x_est = Est_Wind_Height(:,end);
y_est = Est_Wind_Height(:,1);
plot(y_est,x_est, 'o', 'Color', [0 0.5 0.5]);

% Real(Flown) Wind Points
Act_Wind_Height = [W_log.signals.values(N_start+start_filter:N_end,1:2), heights(N_start+start_filter:end)];
Act_Wind_Height = sortrows(Act_Wind_Height, 3);
x_act = Act_Wind_Height(:,end);
y_act = Act_Wind_Height(:,1);
plot(y_act,x_act, 'o', 'Color', [1 0.5 0]);

% Statistical Analysis of the results
% Divide the flown heights in segements, for each clamp the results to the
% lower bound and compute mean and std
ref_heights = linspace(min(heights), max(heights), segments+1);
Est_WindX_meanstd = zeros(size(ref_heights,2)-1,3);
for i = 1:size(ref_heights,2)-1
    h_ind = find(Est_Wind_Height(:,3)>= ref_heights(i) & Est_Wind_Height(:,3)< ref_heights(i+1));
    Est_WindX_meanstd(i,:) = [mean(Est_Wind_Height(h_ind,1)), std(Est_Wind_Height(h_ind,1)), ref_heights(i)];
end
errorbar(Est_WindX_meanstd(:,1),Est_WindX_meanstd(:,3), Est_WindX_meanstd(:,2), 'o-', "horizontal", 'Color', 'black', 'MarkerFaceColor','black', 'LineWidth', 1)

ylim([10 30]), xlabel("Wind Speed [m/s]"), ylabel("Height [m]");
title('Wind Magnitude Profile');
legend('DataPts','Interp', 'Estimated', 'Real','Statistical', 'Location','northwest');