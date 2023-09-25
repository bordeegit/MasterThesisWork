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
Q = diag([1, 1, 1, 5, 5, 5, 30, 30, 30, 5,...  %Q(r,rd,rdd,v)
          5, 1000, 1000, 1000, 1000, 5, 0.1]); %Q(wn,F_L,F_D,cu)
R = diag([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, ...  %R(r,rd)
          0.001, 0.001, 0.001, 1e-15]);                  %R(wr,phir,F_T,delta)

% initial state vector (17)
x = [0.95*Position.signals.values(1,:)'        % r
     1.1*PositionDot.signals.values(1,:)'     % rdot
     [0,0,0]'                             % rdotdot 
     1                                    % Lagrange Multiplier
     [10,10]'                               % W_n 
     [0,0,0]'                                      % F_L vec
     0                                        % F_D
     0];                                        % c_u                                     
us_vec = HLC_input.signals.values;

% Measurements (get index at every loop)

% y = [Position.signals.values(ind,:)'       % r_meas (3)
%      PositionDot.signals.values(ind,:)'    % rdot_meas (3)
%      wr                                    % wr (1)
%      phi_wr                                % phi_wr (1)
%      Forces.signals.values(ind,end)];      % cable traction force |F_T| (1)
%      delta(ind,1)];                        % delta, forcing perp., always 0

f = @(x,us)StateEquations(x, us_vec(1), Ts_10ms, mass, mt_noL);   % State Equations
h = @(x)MeasurementEquations(x, hr, h0);           % Measurement Equations

n = size(x,1);                              % Number of States
P = eye(n);                       % State Covariance
N = 1200;                                 % Number of Steps (12001 max)

xhatV = zeros(N, n);                % Logging estimate  

for ind=1:N
  y = [Position.signals.values(ind,:)'                      % r_meas
       PositionDot.signals.values(ind,:)'                   % rdot_meas
       wr + (-0.5 + rand)                                   % wr with noise
       phi_wr + (-0.1 +(0.2*rand))                          % phi_wr with noise
       Forces.signals.values(ind,end)                       % F_T
       0];                                                  % delta
  [x, P] = StepEKF(f,x,P,h,y,Q,R, diffStepSize, us_vec(ind));   % ekf
  xhatV(ind,:)= x;                                          % Save Estimate
end

