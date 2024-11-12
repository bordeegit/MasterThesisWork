% Flags for Specific plotting
DataFlag = "Kitemill";

% Signals
T_s = 0.01;
pos = [out.logsout.get("x").Values.Data out.logsout.get("y").Values.Data out.logsout.get("h").Values.Data];             % Position in Global frame
posDot = [out.logsout.get("x_dot").Values.Data out.logsout.get("y_dot").Values.Data out.logsout.get("z_dot").Values.Data];    % Speed in Global frame
posDotDot = [0 0 0; diff(posDot)/T_s];                          % Acceleration in Global frame
F_T_norm = out.logsout.get("F_winch").Values.Data;     % Tether Force Magnitude
W = out.logsout.get("wind speed").Values.Data;                  % Absolute Wind (xyz)
Cd_sim = out.logsout.get("CD").Values.Data;                % Cd at each point
Cl_sim = out.logsout.get("CL").Values.Data;                % Cl at each point

% For LS approach  
L_dot = vecnorm(out.logsout.get('xyhdot').Values.Data')'; % (taut) Cable unwinding/winding speed 
beta = reshape(out.logsout.get('alfa_deg').Values.Data, [], 1);  % AoA variation

% Structural Parameters (from functions/parameters.m)
parameters.rho          = 1.225;
parameters.A            = 2.082;
parameters.Cd_l         = 1.2;
parameters.d_l          = 3.5e-3;
parameters.n_l          = 1;
parameters.mk           = 54;
parameters.mt_noL       = 0.74/100; % dyneema-sk99(3.5mm) weight is 0.74kg/100m
parameters.g            = 9.81;
parameters.Cd           = mean(Cd_sim);
parameters.Cl           = mean(Cl_sim);

T_s = 0.01;