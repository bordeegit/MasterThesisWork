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
%   force maybe

diffStepSize = 1e-5;                      % Differentiation Step Size
Q = diag(1, 1, 1, 5, 5, 5, 30, 30, 30,...
         5, 5, 1000, 1000, 1000, 1000, 5, 0.1);
R = diag(0.001, 0.001, 0.001, 0.001, 0.001, 0.001,...
         0.001, 0.001, 0.001, 1e-15);

% initial state vector
x0 = [Position.signals.values(1,:)'        % r
      PositionDot.signals.values(1,:)'     % rdot
      [0,0,0]'                             % rdotdot 
      1                                    % Lagrange Multiplier
      [0,0]'                               % W_n (TODO)
      [0,0,0]'                             % F_L (TODO)
      0                                    % F_D (TODO)
      0];                                  % c_u (TODO)

% Measurements (get index at every loop)
ind = 1;
y = [Position.signals.values(ind,:)'       % r_meas
     PositionDot.signals.values(ind,:)'    % rdot_meas
     wr                                    % wr
     phi_wr                                % phi_wr
     Forces.signals.values(ind,end)        % cable traction force F_T
     delta(ind,1)];                        % delta, forcing perp., always 0

f = @(x)StateEquations(x, diffStepSize);






