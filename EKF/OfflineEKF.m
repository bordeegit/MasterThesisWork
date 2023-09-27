close all 

EKFPrep

%% Summary of elements 

diffStepSize = 1e-5;                      % Differentiation Step Size
Q = diag([1, 1, 1, 5, 5, 5, 30, 30, 30, 5,...  %Q(r,rd,rdd,v)
          5, 1000, 1000, 1000, 1000, 5, 0.1]); %Q(wn,F_L,F_D,cu)
R = diag([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, ...  %R(r,rd)
          0.001, 0.001, 0.001, 1e-15]);                  %R(wr,phir,F_T,delta)

% initial state vector (17) x0
x = [0.9*Position.signals.values(1,:)'        % r
     1.1*PositionDot.signals.values(1,:)'     % rdot
     [1,1,1]'                             % rdotdot 
     0.1                                    % Lagrange Multiplier
     [10,10]'                               % W_n 
     [1,1,1]'                                      % F_L vec
     0.1                                        % F_D
     0.1];                                        % c_u                                     
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
P = 0.1*eye(n);                                 % State Covariance P0
N = 1200;                                 % Number of Steps (12001 max)

xhatV = zeros(N, n);                % Logging estimate  

for ind=1:N
  y = [Position.signals.values(ind,:)'+ [rand,rand,rand]'         % r_meas
       PositionDot.signals.values(ind,:)' + [rand,rand,rand]'     % rdot_meas
       wr + (-0.5 + rand)                                   % wr with noise
       phi_wr + (-0.1 +(0.2*rand))                          % phi_wr with noise
       Forces.signals.values(ind,end) + (-1+2*rand)                        % F_T
       0];                                                  % delta
  [x, P] = StepEKF(f,x,P,h,y,Q,R, diffStepSize, us_vec(ind));   % ekf
  xhatV(ind,:)= x;                                          % Save Estimate
end

