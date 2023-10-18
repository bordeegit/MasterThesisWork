close all 

EKFPrep

% Structural Parameters
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;


% Simulation/Optimization Parameters
parameters.Ts           = Ts_10ms;
W0                      = W_log.signals.values(1,:)';
parameters.r_meas       = Position.signals.values(1,:)';
parameters.rd_meas      = PositionDot.signals.values(1,:)';
parameters.rdd_meas     = PositionDotDot.signals.values(1,:)';
parameters.F_T_norm     = Forces.signals.values(1,end);
parameters.Cl           = Cl.signals.values(1);
parameters.Cd           = Cd.signals.values(1);
parameters.Fl           = Forces.signals.values(1,4:6)'; %!!!
N_opt                   = 1000; % Number of steps to perform optimization
W0_vec                  = [W0';zeros(N_opt-1, 3)];

% Constaints
A = zeros(3,3);
A(3,3) = 1;
b = zeros(3,1);

% Optimization Options
myoptions               = myoptimset;
myoptions.Hessmethod    = 'BFGS';  % Select BFGS or GN
myoptions.gradmethod    = 'CD';
myoptions.graddx        = 2^-17;
myoptions.tolgrad    	= 1e-6;               %default : 1e-6
myoptions.tolfun        = 1e-12;              %default : 1e-12
myoptions.ls_beta       = 0.2;    
myoptions.ls_c          = 0.1;
myoptions.nitermax      = 200;
myoptions.xsequence 	= 'off';
myoptions.display       = 'off';
myoptions.BFGS_gamma 	= 0.1; 
myoptions.GN_sigma      = 0;
myoptions.GN_funF       = @(W)Wind_cost_GN(W,parameters);

codegen Wind_cost -args {W0,parameters} -lang:c++

tic
for i=2:N_opt
    [Wstar,fxstar,k,exitflag,~] = myfmincon(@(W)Wind_cost_mex(W,parameters),W0,A,b,[],[],0,0,myoptions);
    parameters.r_meas       = Position.signals.values(i,:)';
    parameters.rd_meas      = PositionDot.signals.values(i,:)';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
    parameters.F_T_norm     = Forces.signals.values(i,end);
    parameters.Cl           = Cl.signals.values(i);
    parameters.Cd           = Cd.signals.values(i);
    parameters.Fl           = Forces.signals.values(i,4:6)';
    W0_vec(i,:)             = Wstar;
    fprintf("Iteration %d done, Wind is [%f %f %f]\n", i, Wstar(1), Wstar(2), Wstar(3));
end
toc

figure(1)
plot(1:N_opt, W_log.signals.values(1:N_opt,1), 1:N_opt, W0_vec(:,1));
figure(2)
plot(1:N_opt, W_log.signals.values(1:N_opt,2), 1:N_opt, W0_vec(:,2));
figure(3)
plot(1:N_opt, W_log.signals.values(1:N_opt,3), 1:N_opt, W0_vec(:,3));
