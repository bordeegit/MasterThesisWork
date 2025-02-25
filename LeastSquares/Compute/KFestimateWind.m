function [W_est] = KFestimateWind(W_er_measured, l, Ldot, dt, maxStep, turning_mask)
      % Input:
    %   W_er_measured: Wind measurement derived from traction (Nx1 vector)
    %   l: Position unit vector over time (Nx2 matrix)
    %   Ldot: Cable reelout speed (Nx1 vector, can be zeros)
    %   dt: Time step between measurements (scalar, in seconds)
    %   maxStep: final step of the estimation (scalar, in steps)
    %   turning_mask: mask to select invalid points (0 straight, 1 turn)
    %
    % Output:
    %   W_est: Estimated wind matrix over time (Nx2 matrix)
    
    % Ensure inputs are column vectors
    lx = l(:,1); 
    ly = l(:,2);
    Ldot = Ldot(:);
    W_er_measured = W_er_measured(:);
    
     % Handle case where turning_mask is not provided
    if nargin < 7 || isempty(turning_mask)
        turning_mask = zeros(size(W_er_measured)); % Default: no turning
    else
        turning_mask = turning_mask(:); % Ensure column vector
    end
    
    % Number of measurements
    N = maxStep;
    
    % Initialize state and covariance
    x_hat = zeros(2, N);
    x_hat(:,1) = [5; 0];  % Initial guess - adapt based on your knowledge
    
    P = zeros(2, 2, N);
    P(:,:,1) = 5*eye(2);   % Initial uncertainty
    
    % Process noise covariance 
    sigma_w_x = 0.26;  % Wind variation standard deviation [m/s/sqrt(s)]
    sigma_w_y = 0.17;  % Wind variation standard deviation [m/s/sqrt(s)]
    Q = diag([(sigma_w_x^2 * dt),(sigma_w_y^2 * dt)]);  % Scale Q with dt
    
    % Baseline measurement noise (when not turning)
    sigma_meas_base = 0.2;  % Base measurement noise standard deviation [m/s]
    
    % Additional measurement noise during turns
    sigma_meas_turn = 2.0;  % Higher measurement noise during turns [m/s]
    
    % Run Kalman filter
    for k = 2:N
        % Prediction step
        x_minus = x_hat(:,k-1);
        P_minus = P(:,:,k-1) + Q;
        
        % Adaptive measurement noise based on turning status
        if turning_mask(k) == 1
            R_k = sigma_meas_turn^2;
        else
            R_k = sigma_meas_base^2;
        end
        
        % Measurement update step
        H = [lx(k); ly(k)];  % Measurement matrix (2x1 transposed to 1x2)
        z_pred = H' * x_minus - Ldot(k);  % Scalar prediction
        
        % Kalman gain
        K = P_minus * H / (H' * P_minus * H + R_k);  % 2x1 Kalman gain
        
        % Update state and covariance
        x_hat(:,k) = x_minus + K * (W_er_measured(k) - z_pred);
        P(:,:,k) = (eye(2) - K * H') * P_minus;
    end
    
    %KFplotCovariance2

    % Return estimated wind components as column vectors
    W_est = [x_hat(1,:)' x_hat(2,:)' zeros(N,1)];
end