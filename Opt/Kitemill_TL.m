% Flags for Specific plotting
DataFlag = "Kitemill";

% Signals
pos = [out.logsout{8}.Values.Data out.logsout{9}.Values.Data out.logsout{4}.Values.Data];             % Position in Global frame
posDot = out.logsout{35}.Values.Data;       % Speed in Global frame
posDotDot = %TODO; % Acceleration in Global frame
F_T_norm = out.logsout{59}.Values.Data;   % Tether Force Magnitude
W = repmat(out.logsout{16}.Values.Data, size(out.logsout{59}.Values.Data));                  % Absolute Wind (xyz)
Cd_sim = out.logsout{53}.Values.Data;                % Cd at each point
Cl_sim = out.logsout{54}.Values.Data;                % Cl at each point
% Check if Cd and Cl are not inverted 

% Structural Parameters TODO
parameters.rho          = rho;
parameters.A            = area;
parameters.mk           = mass;
parameters.mt_noL       = mt_noL;
parameters.g            = g;
parameters.Cd_l         = CD_Line;
parameters.d_l          = Line_diameter;
parameters.n_l          = n_line;
parameters.Cd           = mean(Cd_sim);
parameters.Cl           = mean(Cl_sim);