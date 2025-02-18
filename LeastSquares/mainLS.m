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

AeroCoeffMode = "mean";      % mean or real 

%% Computation of Equivalent Kite Aerodynamic Efficiency 

if AeroCoeffMode == "mean"
    C_L = mean(Cl_sim)*ones(size(Cl_sim));
    C_D = mean(Cd_sim)*ones(size(Cd_sim)); 
elseif AeroCoeffMode == "real"
    C_L = Cl_sim;
    C_D = Cd_sim;
end

r_l = vecnorm(pos')'; 

% Assuming beta small, cos(beta) ~ 1, so we can neglect it
%   C_Deq = C_D.*(1 + (parameters.n_l*r_l*parameters.d_l*parameters.Cd_l.*cos(beta))./(4*parameters.A*C_D));
C_Deq = C_D.*(1 + (parameters.n_l*r_l*parameters.d_l*parameters.Cd_l)./(4*parameters.A*C_D));

E_eq = C_L./C_Deq; % Equivalent Kite Aerodynamic Efficiency 
%E_eq = 5.3*ones(25001, 1);

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

%% Computation of |W_er| and comparison with real one
% There are 2 methods to compute it (traction or speed method)

F_gr = dot(repmat([0,0,parameters.mk*parameters.g], length(F_T_norm),1),pos./vecnorm(pos,2,2), 2); 
W_er_norm_vec = zeros(length(F_T_norm),2);
C = 0.5*parameters.rho*parameters.A*C_L.*E_eq.^2.*(1+1./E_eq.^2).^(3/2);
%W_er_norm_vec(:,1) = sqrt(abs(F_T_norm./cos(Psi.signals.values) - F_gr - F_app_r)./C);  % Traction approach
W_er_norm_vec(:,1) = sqrt(abs(F_T_norm - F_gr)./C);  % Traction approach
W_er_norm_vec(:,2) = vecnorm(posDot')'./E_eq;  % Speed approach

W_er_norm_real = dot((W - posDot), pos./r_l, 2); %(usare solo vento reale)

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
fprintf("approx_validity: %.3f  %.3f\n",approx_validity);

% Overall, the best choice seems to be the traction


%% Absolute Wind Recovery with Least Squares Approach

% Selection of the type of |W_e,r|
%   1: traction
%   2: speed
typeIndex = 1; 

W_er_norm = W_er_norm_vec(:,typeIndex);  

l = pos./vecnorm(pos,2,2); % kite position versor

blockSize = 4; % Size of the block/window 
initStep = 100; % Starting point, includes a left-truncation
maxStep = 6000; % Ending point (consider T_s=0.02, 6000 steps are 120s)
assert(maxStep <= length(l), "maxStep shouldn't exceed %d", length(l))
noZ = 1;        % Remove the computation of Wz (1 = yes, 0 = no)

%[W_est, iSequence] = blockLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);
tic
[W_est, iSequence] = slidingLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ);
toc

% Smoothing 
W_est_filt = smoothdata(W_est, "movmedian", [500 0], "omitnan"); % Causale


% Performance Factors 
RMSE = rmse(W(iSequence,:), W_est(iSequence,:));
fprintf("RMSE Raw: %f, %f, %f, 2-norm: %f\n", RMSE, norm(RMSE));
RMSE = rmse(W(iSequence,:), W_est_filt(iSequence,:));
fprintf("RMSE Filt: %f, %f, %f, 2-norm: %f\n", RMSE, norm(RMSE));
disp(['For X, Estimated Mean is ' num2str(mean(W_est_filt(:,1), "omitnan"))...
      ' and real mean is ' num2str(mean(W(1:maxStep,1)))])
disp(['For Y, Estimated Mean is ' num2str(mean(W_est_filt(:,2), "omitnan"))...
      ' and real mean is ' num2str(mean(W(1:maxStep,2)))])

load estdata.mat
% Plot Results
timeX = 0:T_s:(maxStep-1)*T_s;
if ~noZ n_subpl = 3; else n_subpl = 2; end
figure, grid on, hold on, %sgtitle("Estimation Results")
subplot(n_subpl,1,1), grid on, hold on,
plot(timeX,W_est(:,1), 'Color', [0.565 0.808 0.98 0.2]);
plot(timeX,W_est_filt(:,1), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(timeX,W(1:maxStep,1),'--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2), %xlim([0 60]), %ylim([-1 15]), hold off
plot(timeX,W0_vec(:,1), 'LineWidth', 2), %xlim([1 6000]), ylim([-1 15]), hold off
legend('Estimated $W_x$', 'Filtered $W_x$', 'Actual $W_x$', 'Opt $W_x$'), ylim([-3 7])
ylabel('Wind Speed (m/s)'), xlabel('Time (s)'), box on
subplot(n_subpl,1,2), grid on, hold on
plot(timeX,W_est(:,2), 'Color', [0.565 0.808 0.98 0.2]);
plot(timeX,W_est_filt(:,2), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(timeX,W(1:maxStep,2),'--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2), %xlim([1 60]), %ylim([-1 5]), hold off
plot(timeX,W0_vec(:,2), 'LineWidth', 2), %xlim([1 6000]), ylim([-1 15]), hold off
legend('Estimated $W_y$', 'Filtered $W_y$', 'Actual $W_y$', 'Opt $W_y$'), ylim([-3 7])
ylabel('Wind Speed (m/s)'), xlabel('Time (s)'), box on
if(~noZ)
subplot(n_subpl,1,3), grid on, hold on,
grid on, hold on, ylim([-5 10])
plot(timeX,W_est(:,3), 'Color', [0.565 0.808 0.98 0.1]);
plot(timeX,W_est_filt(:,3), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(timeX,W(1:maxStep,3),'--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2), %xlim([1 6000]), ylim([-1 15]), hold off
legend('Estimated $W_z$', 'Filtered $W_z$', 'Actual $W_z$'), ylim([-3 7])
ylabel('Wind Speed (m/s)'), xlabel('Time (s)')
end 
linkaxes([subplot(n_subpl,1,1), subplot(n_subpl,1,2)], 'x');  % Link both x and y axes

%exportgraphics(gca,'test2.pdf','ContentType','vector')

% Unfiltered Errors 
error = W(1:maxStep+1,:) - [0,0,0;W_est]; 
error(error(:,1) > 20, 1) = 20; 
error(error(:,1) < -20, 1) = -20; 
error(error(:,2) > 20, 2) = 20;
error(error(:,2) < -20, 2) = -20;  
printTraj(pos, error , initStep, maxStep, "hot")

% Absolute Wind 
W_est_norm = vecnorm(W_est_filt')';
W_norm = vecnorm(W(1:maxStep,:)')';
W0_norm = vecnorm(W0_vec')';
figure,  hold on, grid on
plot(timeX, W_est_norm), plot(timeX, W0_norm), plot(timeX, W_norm, '--'), 
legend('$|W_{est}|$',  '$|W_{OPTest}|$', '$|W|$'), ylabel('Wind Speed (m/s)'), xlabel('Time (s)') 

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