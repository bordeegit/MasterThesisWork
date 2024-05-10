figure; hold on

% Real(Flown) Wind Points
Act_Wind_Height = [W_log.signals.values(N_start:N_end,1:2), heights];
Act_Wind_Height = sortrows(Act_Wind_Height, 3);
x_act = Act_Wind_Height(:,end);
dir_act = atan2d(Act_Wind_Height(:,2),Act_Wind_Height(:,1));
plot(dir_act,x_act, 'o');

Est_Wind_Height = [W0_vec(:,1:2), heights];
Est_Wind_Height = sortrows(Est_Wind_Height, size(Est_Wind_Height,2));
x_est = Est_Wind_Height(:,end);
dir_est = atan2d(Est_Wind_Height(:,2),Est_Wind_Height(:,1));
plot(dir_est,x_est, 'o');