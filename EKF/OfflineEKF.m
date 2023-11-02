close all 

%% Data Preparation 

drag_norm = zeros(12001,1);
for ind = 1:12001
    drag_norm(ind) = norm(Forces.signals.values(ind,7:9));
end


%% EKF Parameters and Execution

diffStepSize = 1e-5;                            % Differentiation Step Size
%Q = diag(0.1*ones(16,1));
Q = diag([1, 1, 1, 5, 5, 5, 30, 30, 30, ... %Q(r,rd,rdd)
          5, 5, 1000, 1000, 1000, 5, 0.1]);        %Q(wn,F_L,F_D,cu)
          
R = diag([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, ...  %R(r,rd)
          0.001, 0.001, 0.001, 1e-15]);                  %R(rdd,delta)
%R = diag(0.0001*ones(10,1));  

% initial state vector x0 (16)
x = [0.5*Position.signals.values(1,:)'          % r
     1.7*PositionDot.signals.values(1,:)'       % rdot
     1.5*PositionDotDot.signals.values(1,:)'    % rdotdot 
     [5,5]'                                     % W_n 
     [1,1,1]'                                   % F_L vec
     1000                                       % F_D
     35];                                       % c_u                                     
us_vec = HLC_input.signals.values;
F_T_vec = Forces.signals.values(:,end);

% True initial state vector, used for P initialization
x0 = [Position.signals.values(1,:)'             % r
      PositionDot.signals.values(1,:)'          % rdot
      PositionDotDot.signals.values(1,:)'       % rdotdot 
      W_log.signals.values(1,1:2)'              % W_n 
      Forces.signals.values(1,4:6)'             % F_L vec
      drag_norm(1)                              % F_D
      37.35];                                     % c_u

% Measurements (get index at every loop)

% y = [Position.signals.values(ind,:)'          % r_meas (3)
%      PositionDot.signals.values(ind,:)'       % rdot_meas (3)
%      PositionDotDot.signals.values(ind,:)'    % rdotdot_meas (3)
%      0];                                      % delta, forcing perp., always 0
    
input = [us_vec(1);F_T_vec(1)];

f = @(x,us)StateEquations(x, input, Ts_10ms, mass, mt_noL);   % State Equations
h = @(x)MeasurementEquations(x);                % Measurement Equations

n = size(x,1);                                  % Number of States
P = diag((x-x0).^2);                            % State Covariance P0
N = 3000;                                      % Number of Steps (12001 max)

xhatV = zeros(N, n);                            % Logging estimate  
zhatV_norm = zeros(N,1);

for ind=1:N
  y = [Position.signals.values(ind,:)';% + 0.5*[rand,rand,rand]';        % r_meas
       PositionDot.signals.values(ind,:)';% + 0.5*[rand,rand,rand]';     % rdot_meas
       PositionDotDot.signals.values(ind,:)';% + 0.5*[rand,rand,rand]';    % rdotdot_meas
       0];                                                               % delta
  [x, P,zhat] = StepEKF(f,x,P,h,y,Q,R, diffStepSize, [us_vec(ind);F_T_vec(ind)]);   % ekf
  xhatV(ind,:)= x;                                % Save Estimate
  zhatV_norm(ind)= norm(zhat);
end

%plot(1:N, zhatV_norm(:))
PlottingEKF