start_filter = 50; % Filter out the first N samples, remove initial transient
segments = 20;     % Number of segments for statistical estimation

f = figure;
labels = [];
f.Position = [400 200 1200 500];

start_filter = max(start_filter - N_start, 0);

%% Plot Wind Profile

subplot(1,2,1);
grid on, hold on

if DataFlag == "SoftKite_CustomLinY"
    % Input Wind Profile, Interpolations
    y_table = 0:0.25:50;
    plot(w_x_height_Data*z_scale,height_Data, 'o', 'Color',[0 0 0.5], 'MarkerFaceColor',[0 0 0.5]);
    x_table = spline(height_Data,w_x_height_Data*z_scale, y_table);
    plot(x_table,y_table, 'Color', [1 0.6 0.2]);
    labels = [labels, {'DataPts'}, {'Interp'}];
end

% Real data, from W = W0/log(Z0/Zr)*log((1/Zr) * Z)
if DataFlag == "SoftKite_LogTrapani"
    hh = 0:0.1:40;
    a = W0_Tr/log(Z0_Tr/Zr_Tr);
    b = 1/Zr_Tr;
    plot(a*log(b*hh), hh)
    labels = [labels, {'Interp'}];
end

% Estimated Wind Points
Est_Wind_Height = [W0_vec(start_filter+1:end,1:2), heights(start_filter+1:end)];
Est_Wind_Height = sortrows(Est_Wind_Height, size(Est_Wind_Height,2));
x_est = Est_Wind_Height(:,end);
y_est = Est_Wind_Height(:,1);
dir_est = atan2d(Est_Wind_Height(:,2),Est_Wind_Height(:,1));
plot(y_est,x_est, 'o', 'Color', [0 0.5 0.5]);
labels = [labels, {'Estimated'}];

% Real(Flown) Wind Points
Act_Wind_Height = [W(N_start+start_filter:N_end,1:2), heights(start_filter+1:end)];
Act_Wind_Height = sortrows(Act_Wind_Height, 3);
x_act = Act_Wind_Height(:,end);
y_act = Act_Wind_Height(:,1);
dir_act = atan2d(Act_Wind_Height(:,2),Act_Wind_Height(:,1));
plot(y_act,x_act, 'o', 'Color', [1 0.5 0]);
labels = [labels, {'Real'}];

% Statistical Analysis of the results
% Divide the flown heights in segements, for each clamp the results to the
% lower bound and compute mean and std
ref_heights = linspace(min(heights), max(heights), segments+1);
Est_WindX_meanstd = zeros(length(ref_heights)-1,3);
Est_WindDir_meanstd = zeros(length(ref_heights)-1,3);
for i = 1:size(ref_heights,2)-1
    h_ind = find(Est_Wind_Height(:,3)>= ref_heights(i) & Est_Wind_Height(:,3)< ref_heights(i+1));
    Est_WindX_meanstd(i,:) = [mean(Est_Wind_Height(h_ind,1)), std(Est_Wind_Height(h_ind,1)), ref_heights(i)];
    Est_WindDir_meanstd(i,:) = [mean(dir_est(h_ind,1)), std(dir_est(h_ind,1)), ref_heights(i)];
end
errorbar(Est_WindX_meanstd(:,1),Est_WindX_meanstd(:,3), Est_WindX_meanstd(:,2), 'o-', "horizontal", 'Color', 'black', 'MarkerFaceColor','black', 'LineWidth', 1)
%ylim([0 40])
xlabel("Wind Speed [m/s]"), ylabel("Height [m]");
labels = [labels, {'Statistical'}];
title('Wind Magnitude Profile');
legend(labels, 'Location','northwest');

% Interpolation of reconstructed points
%xx = 5:0.25:20;
%yy = interp1(Est_WindX_meanstd(:,1), Est_WindX_meanstd(:,3),xx, 'linear');
%plot(xx,yy)

% Using fit (Curve Fitting Toolbox)
% ft = fittype('a*log(b*x)');
% fitted_model = fit(Est_WindX_meanstd(:,1), Est_WindX_meanstd(:,3), ft, 'StartPoint', [1, 1])
% coefficients = coeffvalues(fitted_model);
% a = coefficients(1);
% b = coefficients(2);
% plot(xx,a*log(b*xx))

%% Plot Wind Direction 

subplot(1,2,2);
hold on, grid on

plot(dir_est,x_est, 'o');

plot(dir_act,x_act, 'o');

errorbar(Est_WindDir_meanstd(:,1),Est_WindDir_meanstd(:,3), Est_WindDir_meanstd(:,2), 'o-', "horizontal", 'Color', 'black', 'MarkerFaceColor','black', 'LineWidth', 1);
xlabel("Wind Direction [deg]"), ylabel("Height [m]");
xlim([-10 30])
%ylim([10 30]);
title('Wind Direction Profile');
legend('Estimated', 'Real','Statistical', 'Location','northwest');