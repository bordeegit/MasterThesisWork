clear
close all

set(groot,'DefaultAxesFontSize', 15);
set(groot,'DefaultLineLineWidth', 1.5);
set(groot,'DefaultTextInterpreter', 'Latex');
set(groot,'DefaultAxesTickLabelInterpreter', 'Latex');  
set(groot,'DefaultLegendInterpreter', 'Latex');

%% Load Flight Data & Signals Convertion

% load FlightData\Standard_Step.mat
% SoftKite_TL

% load FlightData\Kitemill_90S.mat
% Kitemill_TL

RealFlightTL

%% Flags 

AeroCoeffMode = "mean";     % mean or real (same in RealFlight) 

ApproxIndex = 1;            % Approx of |W_e,r|    1: traction  2: speed

blockSize = 6;              % Size of the block/window 
initStep = 100;             % Starting point, includes a left-truncation
maxStep = 6000;             % Ending point (consider T_s=0.02, 6000 steps are 120s)

noZ = 1;                    % Remove the computation of Wz (1 = yes, 0 = no)

PrintErrorScaling = false; % Print the plot of the original and scaled error 

%% Computation of Equivalent Kite Aerodynamic Efficiency 

if AeroCoeffMode == "mean"
    C_L = mean(Cl_sim)*ones(size(Cl_sim));
    C_D = mean(Cd_sim)*ones(size(Cd_sim)); 
elseif AeroCoeffMode == "real"
    C_L = Cl_sim;
    C_D = Cd_sim;
end

% Length of the cable, constant in RealFlight
r_l = vecnorm(pos')'; 

% Computation of Equivalent Drag
% Assuming beta small, cos(beta) ~ 1, so we can neglect it
%   C_Deq = C_D.*(1 + (parameters.n_l*r_l*parameters.d_l*parameters.Cd_l.*cos(beta))./(4*parameters.A*C_D));
C_Deq = C_D.*(1 + (parameters.n_l*r_l*parameters.d_l*parameters.Cd_l)./(4*parameters.A*C_D));

E_eq = C_L./C_Deq; 


%% Computation of |W_er| and comparison with real one
% There are 2 methods to compute it (traction or speed method)

F_gr = dot(repmat([0,0,parameters.mk*parameters.g], length(F_T_norm),1),pos./vecnorm(pos,2,2), 2); 
W_er_norm_vec = zeros(length(F_T_norm),2);
C = 0.5*parameters.rho*parameters.A*C_L.*E_eq.^2.*(1+1./E_eq.^2).^(3/2);
%W_er_norm_vec(:,1) = sqrt(abs(F_T_norm./cos(Psi.signals.values) - F_gr - F_app_r)./C);  % Traction approach
W_er_norm_vec(:,1) = sqrt(abs(F_T_norm - F_gr)./C);  % Traction approach
W_er_norm_vec(:,2) = vecnorm(posDot')'./E_eq;  % Speed approach

W_er_norm_real = dot((W - posDot), pos./r_l, 2); %(usare solo vento reale)

% Plot to evaluate goodness of the methods 
figure,
subplot(2,1,1), grid on, hold on
plot(W_er_norm_vec(:,1)), plot(W_er_norm_vec(:,2))
plot(W_er_norm_real)
title('Values'),xlim([0;6000]) 
legend('traction', 'speed', 'real', 'Location','southeast'), hold off
subplot(2,1,2), grid on, hold on
plot(W_er_norm_real - W_er_norm_vec)
title("Differences"), xlim([0;6000])
legend('traction', 'speed', 'Location','southeast' ), hold off
sgtitle('$\mathbf{|\vec{W}_{e,r}|}$') 
approx_validity = vecnorm(W_er_norm_vec - W_er_norm_real);
fprintf("Approx Validity: %.3f  %.3f\n",approx_validity);

% Overall, the best choice seems to be the traction


%% Absolute Wind Recovery with Least Squares Approach, Different Estimations

% Selection of the approximation
W_er_norm = W_er_norm_vec(:,ApproxIndex);  

% kite position versor
l = pos./vecnorm(pos,2,2); 
assert(maxStep <= length(l), "maxStep shouldn't exceed %d", length(l))

% Basic Estimation with LS 
%[W_est, iSequence] = blockLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);
tic
W_est = slidingLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);
toc

% a) Smoothing of Basic Estimation with moving median (casual)
filtWindow = [500 0]; % Filter window [b f], b element before, f after
W_est_filt = smoothdata(W_est, "movmedian", filtWindow, "omitnan"); 

% b) Filtering of Basic Estimation considering only crosswind flight (CWF)
CWF_TH_UPP = 0.5;
CWF_TH_LOW = -0.5;
[W_est_CWF, CWF_mask] = CWF_Selection(W_est, phi_dot, CWF_TH_UPP, CWF_TH_LOW);

% c) CWF Estimation, sliding only on blocks where all elements are valid 
W_CWF = slidingLS_CWF(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ,CWF_mask);


% d) Estimation using Optimization appraoch (W0_vec)
load estdata.mat
W_OPT = W0_vec;

% e) Real wind measurement
W_real = W(2:maxStep+1,:);

clearvars W W0_vec
W_matrices = {W_est_filt, W_est_CWF, W_CWF, W_OPT, W_real}; 
Names_matrices = {'Est_filt', strcat('Est_CWF_',num2str(CWF_TH_UPP)), 'CWF', 'OPT', 'Real'};


%% Performance Evalation or Plot results 

% Performance Factors wrt real wind (note: at ground)
PerformanceEvaluation(W_matrices, Names_matrices);


% Plot Results TODO: move this in a function, maybe loop over W_matrices
timeX = 0:T_s:(maxStep-1)*T_s;
if ~noZ n_subpl = 3; else n_subpl = 2; end
figure, grid on, hold on, %sgtitle("Estimation Results")
subplot(n_subpl,1,1), grid on, hold on,
plot(timeX,W_est(:,1), 'Color', [0.565 0.808 0.98 0.2]);
plot(timeX,W_est_filt(:,1), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(timeX,W_est_CWF(:,1), 'LineWidth', 2);
plot(timeX,W_CWF(:,1), 'LineWidth', 2);
plot(timeX,W_real(:,1),'--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2), %xlim([0 60]), %ylim([-1 15]), hold off
plot(timeX,W_OPT(:,1), 'LineWidth', 2), %xlim([1 6000]), ylim([-1 15]), hold off
legend('Estimated $W_x$', 'Filtered $W_x$', 'Crosswind $W_x$', 'Crosswind Full $W_x$', 'Actual $W_x$', 'Opt $W_x$'), ylim([-3 7])
ylabel('Wind Speed (m/s)'), xlabel('Time (s)'), box on
subplot(n_subpl,1,2), grid on, hold on
plot(timeX,W_est(:,2), 'Color', [0.565 0.808 0.98 0.2]);
plot(timeX,W_est_filt(:,2), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(timeX,W_est_CWF(:,2), 'LineWidth', 2);
plot(timeX,W_CWF(:,2), 'LineWidth', 2);
plot(timeX,W_real(:,2),'--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2), %xlim([1 60]), %ylim([-1 5]), hold off
plot(timeX,W_OPT(:,2), 'LineWidth', 2), %xlim([1 6000]), ylim([-1 15]), hold off
legend('Estimated $W_y$', 'Filtered $W_y$', 'Crosswind $W_y$', 'Crosswind Full $W_x$', 'Actual $W_y$', 'Opt $W_y$'), ylim([-3 7])
ylabel('Wind Speed (m/s)'), xlabel('Time (s)'), box on
if(~noZ)
subplot(n_subpl,1,3), grid on, hold on,
grid on, hold on, ylim([-5 10])
plot(timeX,W_est(:,3), 'Color', [0.565 0.808 0.98 0.1]);
plot(timeX,W_est_filt(:,3), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(timeX,W_real(:,3),'--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2), %xlim([1 6000]), ylim([-1 15]), hold off
legend('Estimated $W_z$', 'Filtered $W_z$', 'Actual $W_z$'), ylim([-3 7])
ylabel('Wind Speed (m/s)'), xlabel('Time (s)')
end 
linkaxes([subplot(n_subpl,1,1), subplot(n_subpl,1,2)], 'x');  % Link both x and y axes

%exportgraphics(gca,'test2.pdf','ContentType','vector')

% Visualization of zones where spikes accour 
error = W_real - W_est;
error_scaled = errorScaling(error, 'log', PrintErrorScaling);
printTraj(pos, error_scaled , initStep, maxStep, "hot")

% Absolute Wind magnitude
W_est_filt_norm = vecnorm(W_est_filt')';
W_real_norm = vecnorm(W_real')';
W_CWF_norm = vecnorm(W_CWF')';
W_OPT_norm = vecnorm(W_OPT')';
figure,  hold on, grid on
plot(timeX, W_est_filt_norm), plot(timeX, W_OPT_norm), plot(timeX, W_CWF_norm), plot(timeX, W_real_norm, '--'), 
legend('$|W_{filt}|$',  '$|W_{OPT}|$', '$|W_{CWF}|$', '$|W_{real}|$'), ylabel('Wind Speed (m/s)'), xlabel('Time (s)')

% Absolute Wind direction  
dir_W_est_filt = atan2d(W_est_filt(:,2),W_est_filt(:,1));
dir_W = atan2d(W_real(:,2),W_real(:,1));
dir_W_CWF_full = atan2d(W_CWF(:,2),W_CWF(:,1));
dir_W_OPT = atan2d(W_OPT(:,2),W_OPT(:,1));
figure,  hold on, grid on
plot(timeX, dir_W_est_filt), plot(timeX, dir_W_OPT), plot(timeX, dir_W_CWF_full), plot(timeX, dir_W, '--'), 
legend('$\angle W_{filt}$',  '$\angle W_{OPT}$', '$\angle W_{CWF}$', '$\angle W_{real}$'), ylabel('Wind Direction (deg)'), xlabel('Time (s)')

%% Hyperparameter optimization?

%% Alternative approach and other code
% The method works well for constant wind. Then, another idea can be 
% to divide the flight data into heights segmentes, where we can suppose
% constant wind, and estimate there. 
% We can then reconstruct the wind profile by merging the estimation in
% the various heights
%
% Block size wrt loop Period
% loop_period = diff(find(diff(sign(diff(pos(initStep:maxStep,2))))==-2));
% figure, grid on, hold on,
% plot(loop_period), yline(blockSize, 'r--', 'LineWidth', 4), hold off

%% Notes on the approach with beta
% NOTE: The introduction of beta gives minimal difference, the 
%       difference is 2 orders of magnitude less than the absolute value 
%       However, beta is used in the second method
%
% Alternative and equivalent computation of beta (from definition)
%   This ofc cannot be used in practice, since we don't have We (bc of W)
%We = W - posDot;
%beta = pi/2 - acos(dot(We, pos, 2)./(vecnorm(We')'.*vecnorm(pos')'));
%
% This cannot be used bc we cannot have beta, we would need to know the AoA
% which is not know 
% E_eq_2 = cos(beta)./sin(beta); % Equivalent to 1./tan(beta)
% 
% figure, 
% subplot(2,1,1), grid on, hold on;
% plot(E_eq_1), plot(E_eq_2), xlim([0;6000])
% legend('base method','approx method', 'Location','southeast' ), hold off
% subplot(2,1,2), 
% plot(E_eq_1-E_eq_2), grid on, xlim([0;6000]);
% legend('Difference b/w methods')
% sgtitle('Equivalent Kite Aerodyn Efficiency')

%% Computation of F_app_r and testing

% Checked that F_gr is the same as the a priori projection below
%F_grav_r = parameters.mk*parameters.g.*cos(th);

% Apparent force projected on the cable
%F_app_r = parameters.mk*(r_l.*thd.^2 + r_l.*phid.^2.*sin(th).^2);

% Both forces F_gr and F_app_r are negligible wrt to F_T_norm, 1000 vs 10