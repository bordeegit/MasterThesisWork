close all 

EKFPrep

%% Summary of elements 
% Position.signals.values(:,:) position r
% PositionDot.signals.values(:,:) velocity rdot 
% wr wind speed at reference height zr
% phi_wr direction of wind at reference height zr
% Forces.signals.values(:,end) cable traction force F_T
% delta is zero always (see paper)

% F_L and F_D seems hard to compute, especially if we are in global coords
%   We can just skip them for now (differently from the paper), since we 
%   don't care about them (also c_u will be skipped) just estimate total
%   force maybe (this is not correct since we have apparent forces, but we
%   can try)

diffStepSize = 1e-5;                      % Differentiation Step Size
Q = diag([1, 1, 1, 5, 5, 5, 30, 30, 30,...
         5, 5, 1000, 1000, 1000, 1000]);
R = diag([0.001, 0.001, 0.001, 0.001, 0.001, 0.001,...
         0.001, 0.001, 0.001]);

% initial state vector
x = [Position.signals.values(1,:)'        % r
      PositionDot.signals.values(1,:)'     % rdot
      [0,0,0]'                             % rdotdot 
      1                                    % Lagrange Multiplier
      [10,10]'                               % W_n (TODO)
      Forces.signals.values(1,1:3)'+Forces.signals.values(1,4)']; % F_LDg (Maybe wrong ref sys?)

% Measurements (get index at every loop)
% NOTE: since we don't compute F_L directly, we cannot use delta (for now)
% y = [Position.signals.values(ind,:)'       % r_meas
%      PositionDot.signals.values(ind,:)'    % rdot_meas
%      wr                                    % wr
%      phi_wr                                % phi_wr
%      Forces.signals.values(ind,end)];      % cable traction force F_T
%      %delta(ind,1)];                        % delta, forcing perp., always 0

f = @(x)StateEquations(x, Ts_10ms, mass);   % State Equations
h = @(x)MeasurementEquations(x, hr, h0);           % Measurement Equations

n = size(x,1);                              % Number of States
P = eye(n);                       % State Covariance
N = 12001;                                 % Number of Steps

xhatV = zeros(N, n);                % Logging estimate  

for ind=1:N
  y = [Position.signals.values(ind,:)'          % r_meas
       PositionDot.signals.values(ind,:)'       % rdot_meas
       wr                                       % wr
       phi_wr                                   % phi_wr
       Forces.signals.values(ind,end)];         % F_T
  [x, P] = StepEKF(f,x,P,h,y,Q,R, diffStepSize);% ekf
  xhatV(ind,:)= x;                              % Save Estimate
end

