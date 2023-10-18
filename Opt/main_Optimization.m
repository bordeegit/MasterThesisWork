close all 

EKFPrep

% Structural Parameters
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;


% Simulation/Optimization Parameters
parameters.Ts           = Ts_10ms; % Not currently used
z0                      = [W_log.signals.values(1,:)'
                           1;1;1];
parameters.r_meas       = Position.signals.values(1,:)';
parameters.rd_meas      = PositionDot.signals.values(1,:)';
parameters.rdd_meas     = PositionDotDot.signals.values(1,:)';
parameters.F_T_norm     = Forces.signals.values(1,end);
parameters.Cl           = Cl.signals.values(1);
parameters.Cd           = Cd.signals.values(1);
%parameters.Fl           = Forces.signals.values(1,4:6)'; %!!!
N_opt                   = 1000; % Number of steps to perform optimization
W0_vec                  = [W0';zeros(N_opt-1, 3)];
parameters.Q            = diag(100*ones(3,1));

% Constaints
A = zeros(6,6);
A(3,3) = 1;
b = zeros(6,1);

% Optimization Options
myoptions               = myoptimset;
myoptions.Hessmethod    = 'GN';  % Select BFGS or GN
myoptions.gradmethod    = 'CD';
myoptions.graddx        = 2^-17;
myoptions.tolgrad    	= 1e-6;               %default : 1e-6
myoptions.tolfun        = 1e-12;              %default : 1e-12
myoptions.ls_beta       = 0.2;    
myoptions.ls_c          = 0.1;
myoptions.nitermax      = 200;
myoptions.xsequence 	= 'off';
myoptions.display       = 'off';  % or Iter
myoptions.BFGS_gamma 	= 0.1; 
myoptions.GN_sigma      = 0;
%myoptions.GN_funF       = @(z)Wind_cost_GN(z,parameters);

%codegen Wind_cost -args {z0,parameters} -lang:c++
codegen Wind_cost_GN -args {z0,parameters} -lang:c++

tic
for i=2:N_opt
    [zstar,fxstar,k,exitflag,~] = myfmincon(@(z)Wind_cost_GN_mex(z,parameters),z0,A,b,[],[],1,0,myoptions);
    %z0                      = zstar;
    parameters.r_meas       = Position.signals.values(i,:)';
    parameters.rd_meas      = PositionDot.signals.values(i,:)';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
    parameters.F_T_norm     = Forces.signals.values(i,end);
    parameters.Cl           = Cl.signals.values(i);
    parameters.Cd           = Cd.signals.values(i);
    %parameters.Fl           = Forces.signals.values(i,4:6)';
    W0_vec(i,:)             = zstar(1:3);
    fprintf("Iteration %d done, Wind is [%f %f %f], Zl is [%f %f %f]\n", i, zstar(1), zstar(2), zstar(3),zstar(4), zstar(5), zstar(6));
end
toc

figure(1)
plot(1:N_opt, W_log.signals.values(1:N_opt,1), 1:N_opt, W0_vec(:,1));
figure(2)
plot(1:N_opt, W_log.signals.values(1:N_opt,2), 1:N_opt, W0_vec(:,2));
figure(3)
plot(1:N_opt, W_log.signals.values(1:N_opt,3), 1:N_opt, W0_vec(:,3));
