% Flags for Specific plotting
DataFlag = "SoftKite";

% Signals
T_s = 0.01;
pos = Position.signals.values;             % Position in Global frame
posDot = PositionDot.signals.values;       % Speed in Global frame
posDotDot = PositionDotDot.signals.values; % Acceleration in Global frame
F_T_norm = Forces.signals.values(:,end);   % Tether Force Magnitude
W = W_log.signals.values;                  % Absolute Wind (xyz)
Cd_sim = Cd.signals.values;                % Cd at each point
Cl_sim = Cl.signals.values;                % Cl at each point

% For LS approach
L_dot = statesdot.signals.values(:,5);     % Cable unwinding/winding speed
%beta = alpha.signals.values - alpha_0;     % AoA variation
th = pi/2 - theta.signals.values;
thd = - thetadot.signals.values;
phi_dot = phidot.signals.values;

% Structural Parameters
parameters.rho          = rho;
parameters.A            = area;
parameters.Cd_l         = CD_Line;
parameters.d_l          = Line_diameter;
parameters.n_l          = n_line;
parameters.mk           = mass;
parameters.mt_noL       = mt_noL;
parameters.g            = g;
parameters.Cd           = mean(Cd_sim);
parameters.Cl           = mean(Cl_sim);