clear
close all

set(groot,'DefaultAxesFontSize', 15);
set(groot,'DefaultLineLineWidth', 1.5);
set(groot,'DefaultTextInterpreter', 'Latex');
set(groot,'DefaultAxesTickLabelInterpreter', 'Latex');  
set(groot,'DefaultLegendInterpreter', 'Latex');

%% Load Flight Data & Signals Convertion
%load FlightData\Standard_Step.mat
%SoftKite_TL

% Additionally, in TL, get 
%   - L_dot
%   - beta 

load FlightData\Kitemill_Full_180S.mat
Kitemill_TL

%% Flags 

AeroCoeffMode = "mean";      % mean or real 

%% Computation of Equivalent Kite Aerodynamic Efficiency 

% There are 2 ways to compute beta (AoA variation)
% There are 2 ways to compute E_eq (base and approx)
if AeroCoeffMode == "mean"
    C_L = mean(Cl_sim)*ones(size(Cl_sim));
    C_D = mean(Cd_sim)*ones(size(Cd_sim)); 
elseif AeroCoeffMode == "real"
    C_L = Cl_sim;
    C_D = Cd_sim;
end

%n_line = 1;
r_l = vecnorm(pos')'; % equivalent to states.signals.values(:,5)

% Alternative and equivalent computation of beta (from definition)
%   This ofc cannot be used in practice, since we don't have We (bc of W)
%We = W - posDot;
%beta = pi/2 - acos(dot(We, pos, 2)./(vecnorm(We')'.*vecnorm(pos')'));

C_Deq = C_D.*(1 + (parameters.n_l*r_l*parameters.d_l*parameters.Cd_l.*cos(beta))./(4*parameters.A*C_D));

% NOTE: The introduction of beta gives minimal difference, the 
%       difference is 2 orders of magnitude less than the absolute value 
%       However, beta is used in the second method

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

%% Computation of |W_er| and comparison with real one
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
plot(W_er_norm_vec - W_er_norm_real)
% TODO: Ha senso fare grafico delle 4 differenze tra approx e reali?

%plot(W_er_norm_vec(:,1) - W_er_norm_vec(:,2))
%plot(W_er_norm_vec(:,3) - W_er_norm_vec(:,4))
%plot(W_er_norm_vec(:,1) - W_er_norm_vec(:,3))
%plot(W_er_norm_vec(:,2) - W_er_norm_vec(:,4))
title("Differences"), xlim([0;6000])
legend('base E eq','approx E eq', 'traction', 'speed', 'Location','southeast' ), hold off
sgtitle('$\mathbf{|\vec{W}_{e,r}|}$', 'Interpreter','latex') 
approx_validity = vecnorm(W_er_norm_vec - W_er_norm_real);
fprintf("approx_validity: %.3f  %.3f  %.3f  %.3f\n",approx_validity);

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

typeIndex = 2; 
% Selection of the type of |W_e,r|
%   1: base + traction
%   2: base + speed
%   3: approx + traction
%   4: approx + speed

W_er_norm = W_er_norm_vec(:,typeIndex);  

l = pos./vecnorm(pos,2,2); % kite position versor

blockSize =1700; % Size of the block/window 
initStep = 500; % Starting point, includes a left-truncation
maxStep = 10000; % Shouldn't exceed length(l)
noZ = 0;        % Remove the computation of Wz (1 = yes, 0 = no)

%[W_est, iSequence] = blockLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);

[W_est, iSequence] = slidingLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);

% Performance Factors 
RMSE = rmse(W(iSequence,:), W_est(iSequence,:));
fprintf("RMSE: %f, %f, %f, 2-norm: %f\n", RMSE, norm(RMSE));

% Block size wrt loop Period
loop_period = diff(find(diff(sign(diff(pos(initStep:maxStep,2))))==-2));
figure, grid on, hold on,
plot(loop_period), yline(blockSize, 'r--', 'LineWidth', 4), hold off

% Results
figure, grid on, hold on,
subplot(3,1,1), grid on, hold on
plot(W_est(:,1));
plot(W(1:maxStep,1),'--'), %xlim([1 6000]), ylim([-1 15]), hold off
subplot(3,1,2), grid on, hold on
plot(W_est(:,2));
plot(W(1:maxStep,2),'--'), %xlim([1 6000]), ylim([-1 5]), hold off
subplot(3,1,3), grid on, hold on
plot(W_est(:,3));
plot(W(1:maxStep,3),'--'), %xlim([1 6000]), hold off
sgtitle("Estimation Results")

%% Filtering 
% Initial testing with simple threshold approach
%   spike happen around the real/correct value, so a simple 
%   |W_est_vals| < th can't work
%
% Used Matlab's DataCleaner App on VALUES, trial&error for params
%   - filloutliers: remove spikes (moving median) and fill (lin interp) 
%   - smoothdata: smooth data (moving mean with moving window)
% W_est_trunc = W_est_vals(iFilterStart:end,:); % Skip first second 
% W_est_table = array2table(W_est_trunc, 'VariableNames',{'Wx','Wy','Wz'});
% W_est_table = filloutliers(W_est_table,"linear","movmedian",20,"DataVariables",["Wx","Wy"]);
% W_est_table = smoothdata(W_est_table,"movmean",20,"DataVariables",["Wx","Wy"]);
% W_est_vals_fil = table2array(W_est_table);
% 
% figure;
% subplot(3,1,1)
% plot(W_est_trunc(:,1), 'b-'), hold on, grid on,
% plot(W_est_vals_fil(:,1), 'r-');
% legend('Original','Filtered'), hold off;
% subplot(3,1,2)
% plot(W_est_trunc(:,2), 'b-'), hold on, grid on,
% plot(W_est_vals_fil(:,2), 'r-');
% legend('Original','Filtered'), hold off;
% subplot(3,1,3)
% plot(W_est_trunc(:,3), 'b-'), hold on, grid on,
% plot(W_est_vals_fil(:,3), 'r-');
% legend('Original','Filtered'), hold off;
% sgtitle("Filtered values results")
% 
% W_est_rec_fil = kron(W_est_vals_fil, ones(blockSize,1));
% 
% figure, 
% subplot(2,1,1), grid on,hold on
% plot(W_est_rec(N_filter:end,1),'--', 'LineWidth', 1);
% plot(W_est_rec_fil(:,1));
% plot(W(N_filter:end,1)),
% legend('Estimated', 'Filtered', 'Real'), ylim([-5;20]), hold off
% subplot(2,1,2),grid on, hold on
% plot(W_est_rec(N_filter:end,2), '--', 'LineWidth', 1);
% plot(W_est_rec_fil(:,2));
% plot(W(N_filter:end,2)),
% legend('Estimated', 'Filtered', 'Real'), ylim([-5;20]), hold off
% sgtitle("Filtered reconstructed results")
% 
% meanX_fil = mean(W_est_rec_fil(:,1))
% meanY_fil = mean(W_est_rec_fil(:,2))
% meanZ_fil = mean(W_est_rec_fil(:,3))
%% Hyperparameter optimization?

%% Alternative approach
% The method works well for constant wind. Then, another idea can be 
% to divide the flight data into heights segmentes, where we can suppose
% constant wind, and estimate there. 
% We can then reconstruct the wind profile by merging the estimation in
% the various heights
