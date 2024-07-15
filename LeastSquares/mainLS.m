clear
close all

set(groot,'DefaultAxesFontSize', 15);
set(groot,'DefaultLineLineWidth', 1.5);
set(groot,'DefaultTextInterpreter', 'Latex');
set(groot,'DefaultAxesTickLabelInterpreter', 'Latex');  
set(groot,'DefaultLegendInterpreter', 'Latex');

%% Load Flight Data & Signals Convertion
load FlightData\Standard_LinY.mat
SoftKite_TL

%% Computation of Equivalent Kite Aerodynamic Efficiency 
% There are 2 ways to compute beta (AoA variation)
% There are 2 ways to compute E_eq (base and approx)
C_L = Cl_sim;
C_D = Cd_sim; 
r_l = vecnorm(pos')'; % equivalent to states.signals.values(:,5)
beta = alpha.signals.values - alpha_0;
C_Deq = C_D.*(1 + (n_line*r_l*parameters.d_l*parameters.Cd_l.*cos(beta))./(4*parameters.A*C_D));

% Alternative and equivalent computation of beta (from definition)
%   We = W - posDot;
%   beta_alt = pi/2 - acos(dot(We, pos, 2)./(vecnorm(We')'.*vecnorm(pos')'));

% NOTE: The introduction of beta gives minimal difference, the 
%       difference is 2 orders of magnitude less than the absolute value 

E_eq_1 = C_L./C_Deq;
E_eq_2 = cos(beta)./sin(beta); % Equivalent to 1./tan(beta)

figure, 
subplot(2,1,1), grid on, hold on;
plot(E_eq_1), plot(E_eq_2), xlim([0;6000])
legend('base method','approx method', 'Location','southeast' ), hold off
subplot(2,1,2), 
plot(E_eq_1-E_eq_2), grid on, xlim([0;6000]);
legend('Difference b/w methods')
sgtitle('Equivalent Kite Aerodyn Efficiency')

% The 2 methods seems to match up well, the base should be less prone to
% wrong assumptions, and should be more constant 

%% Computation of |W_er| 
% There are 2 methods to compute it (traction or speed method)

W_er_norm_vec = zeros(length(F_T_norm),4);
i = 1;
for E_eq = [E_eq_1,E_eq_2]
    C = 0.5*parameters.rho*parameters.A*C_L.*E_eq.^2.*(1+1./E_eq.^2).^(3/2);
    W_er_norm_vec(:,i) = sqrt(F_T_norm./C); % Traction approach
    W_er_norm_vec(:,i+1) = vecnorm(posDot')'./E_eq;  % Speed approach
    i = i+2;
end

W_er_norm_real = dot((W - posDot), pos./vecnorm(pos,2,2), 2);

figure
subplot(2,1,1), grid on, hold on
plot(W_er_norm_vec(:,1)), plot(W_er_norm_vec(:,2))
plot(W_er_norm_vec(:,3)), plot(W_er_norm_vec(:,4))
plot(W_er_norm_real)
title('Values'),xlim([0;6000]) 
legend('traction base', 'speed base', 'traction approx', 'speed approx', 'real', 'Location','southeast'), hold off
subplot(2,1,2), grid on, hold on
plot(W_er_norm_vec(:,1) - W_er_norm_vec(:,2))
plot(W_er_norm_vec(:,3) - W_er_norm_vec(:,4))
plot(W_er_norm_vec(:,1) - W_er_norm_vec(:,3))
plot(W_er_norm_vec(:,2) - W_er_norm_vec(:,4))
title("Differences"), xlim([0;6000])
legend('base E eq','approx E eq', 'traction', 'speed', 'Location','southeast' ), hold off
sgtitle('$\mathbf{|\vec{W}_{e,r}|}$', 'Interpreter','latex')

% The computation of W_e,r with the unapproximaed method of computation of 
% E_eq produces a lot less spikes, as it can be deduced from the plots of 
% E_eq_1 and E_eq_2 

% Considering the same Equivalent aerodynamic efficiency, the 2 methods of 
% computation of W_e,r (traction and speed) produce very different results;
% in particular, the speed approach more spikes, and in general has a more 
% erratic behivour

% The best one can be found with min(vecnorm(W_er_norm_vec - W_er_norm_real))
% Overall, the best choice seems to be approx E_eq + traction


%% Absolute Wind Recovery with Least Squares Approach

typeIndex = 3; 
% Selection of the type of |W_e,r|
%   1: base + traction
%   2: base + speed
%   3: approx + traction
%   4: approx + speed

W_er_norm = W_er_norm_vec(:,typeIndex);  
L_dot = statesdot.signals.values(:,5); % unwinding/winding speed
l = pos./vecnorm(pos,2,2); % kite position versor

blockSize = 40;
initStep = 200;
maxStep = 6001; % Shouldn't exceed length(l)
noZ = 1;        % Remove the computation of Wz (1 = yes, 0 = no)

%[W_est_rec, W_est_vals, iSequence] = blockLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);

[W_est_rec, W_est_vals, iSequence] = slidingLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);


figure,
subplot(3,1,1),
plot(W_est_vals(:,1)),  grid on;
subplot(3,1,2),
plot(W_est_vals(:,2)),  grid on;
subplot(3,1,3)
plot(W_est_vals(:,3)),  grid on;
sgtitle("Unfiltered Estimation Values")

%% Left-Trucation 

% Due to very high spikes at the start of the simulation, I decided to
% left-truncate the results

% N_filter = 200; % Steps to skip in the comparison (100 steps = 1 sec)
% 
% iFilterStart = 1+floor(N_filter/blockSize);
% iSequence = iSequence(1,iFilterStart:end);
% figure,
% subplot(2,1,1), grid on, hold on
% plot(W_est_rec(N_filter:end,1));
% plot(W(N_filter:end,1),'--'), hold off
% subplot(2,1,2), grid on, hold on
% plot(W_est_rec(N_filter:end,2));
% plot(W(N_filter:end,2),'--'), hold off
% sgtitle("Left-Truncated Estimation Vector")
% 
% % Performance Factors 
% diffValues = abs(norm(W_est_vals(iFilterStart:end,1) - W(iSequence,1)))
% diffRec    = abs(norm(W_est_rec(N_filter:5998,1)-W(N_filter:5998,1)))
% meanX      = mean(W_est_rec(100:end,1))
% meanY      = mean(W_est_rec(100:end,2))
% meanZ      = mean(W_est_rec(100:end,3))

figure, grid on, hold on,
subplot(2,1,1), grid on, hold on
plot(W_est_rec(:,1));
plot(W(1:maxStep,1),'--'), hold off
subplot(2,1,2), grid on, hold on
plot(W_est_rec(:,2));
plot(W(1:maxStep,2),'--'), hold off
sgtitle("Left-Truncated Estimation Vector")

stop 
%% Filtering 

% Initial testing with simple threshold approach
%   spike happen around the real/correct value, so a simple 
%   |W_est_vals| < th can't work

% Used Matlab's DataCleaner App on VALUES, trial&error for params
%   - filloutliers: remove spikes (moving median) and fill (lin interp) 
%   - smoothdata: smooth data (moving mean with moving window)
W_est_trunc = W_est_vals(iFilterStart:end,:); % Skip first second 
W_est_table = array2table(W_est_trunc, 'VariableNames',{'Wx','Wy','Wz'});
W_est_table = filloutliers(W_est_table,"linear","movmedian",20,"DataVariables",["Wx","Wy"]);
W_est_table = smoothdata(W_est_table,"movmean",20,"DataVariables",["Wx","Wy"]);
W_est_vals_fil = table2array(W_est_table);

figure;
subplot(3,1,1)
plot(W_est_trunc(:,1), 'b-'), hold on, grid on,
plot(W_est_vals_fil(:,1), 'r-');
legend('Original','Filtered'), hold off;
subplot(3,1,2)
plot(W_est_trunc(:,2), 'b-'), hold on, grid on,
plot(W_est_vals_fil(:,2), 'r-');
legend('Original','Filtered'), hold off;
subplot(3,1,3)
plot(W_est_trunc(:,3), 'b-'), hold on, grid on,
plot(W_est_vals_fil(:,3), 'r-');
legend('Original','Filtered'), hold off;
sgtitle("Filtered values results")

W_est_rec_fil = kron(W_est_vals_fil, ones(blockSize,1));

figure, 
subplot(2,1,1), grid on,hold on
plot(W_est_rec(N_filter:end,1),'--', 'LineWidth', 1);
plot(W_est_rec_fil(:,1));
plot(W(N_filter:end,1)),
legend('Estimated', 'Filtered', 'Real'), ylim([-5;20]), hold off
subplot(2,1,2),grid on, hold on
plot(W_est_rec(N_filter:end,2), '--', 'LineWidth', 1);
plot(W_est_rec_fil(:,2));
plot(W(N_filter:end,2)),
legend('Estimated', 'Filtered', 'Real'), ylim([-5;20]), hold off
sgtitle("Filtered reconstructed results")

meanX_fil = mean(W_est_rec_fil(:,1))
meanY_fil = mean(W_est_rec_fil(:,2))
meanZ_fil = mean(W_est_rec_fil(:,3))