close all 

% Structural Parameters
parameters.rho = rho;
parameters.A = area;
parameters.mk = mass;
parameters.mt_noL = mt_noL;
parameters.g = g;


% Simulation/Optimization Parameters
parameters.Ts           = Ts_10ms; % Not currently used
z0                      = [W_log.signals.values(1,:)';
                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
%z0                      = [8;3;0; %W_log.signals.values(1,:)';
%                           0.9; 0.1; sqrt(1-0.9^2-0.1^2)]; %1/sqrt(3)*ones(3,1)];
Nop                     = size(z0,1);
parameters.r_meas       = Position.signals.values(1,:)';
parameters.rd_meas      = PositionDot.signals.values(1,:)';
parameters.rdd_meas     = PositionDotDot.signals.values(1,:)';
parameters.F_T_norm     = Forces.signals.values(1,end);
parameters.Cl           = Cl.signals.values(1);
parameters.Cd           = Cd.signals.values(1);
N_opt                   = 500; % Number of steps to perform optimization
W0_vec                  = [z0(1:3)';zeros(N_opt-1, 3)];
parameters.Q            = 1e6*diag(ones(3,1));
parameters.Qw           = 1e3*diag(ones(6,1));
parameters.zold         = z0;
parameters.deltaNorm    = 0.01;   % Delta for Zl nonlin. ineq. constrain
parameters.deltaOrth    = 0.1;   % Delta for Zl*wa nonlin. ineq. constrain
parameters.deltaWind    = [5;1]; % Delta for diff btw Zl_x, Zl_y and old

% Equality Constaints (W_z = 0)
A = zeros(Nop,Nop);
A(3,3) = 1;
b = zeros(Nop,1);

% Inequality Constrains (bounds on W_x, W_y)
C = zeros(2*Nop,Nop);
C(1,1) = -1;
C(2,2) = -1;
C(7,1) = 1;
C(8,2) = 1;
d = [-15;-10; zeros(4,1);
     5;-2; zeros(4,1)];

% Numbers of nonlinera constraints
p = 0;
q = 4+4;

% Optimization Options
myoptions               = myoptimset;
myoptions.Hessmethod    = 'GN';  % Select BFGS or GN
myoptions.gradmethod    = 'CD';
myoptions.graddx        = 2^-17;              %default (CD) :  2^-17
myoptions.tolgrad    	= 1e-6;               %default : 1e-6
myoptions.tolfun        = 1e-12;              %default : 1e-12  
myoptions.ls_beta       = 0.2;%0.2;                %default : 0.8       
myoptions.ls_c          = 0.1;                %default : 0.1
myoptions.ls_nitermax   = 20; %100                %default : 20
myoptions.nitermax      = 100; %200;                %default : 50
myoptions.xsequence 	= 'off';
myoptions.display       = 'off';  % or Iter
myoptions.BFGS_gamma 	= 0.1;                %default : 1e-1
myoptions.GN_sigma      = 0;
myoptions.GN_funF       = @(z)Wind_cost_GN(z,parameters);

%codegen Wind_cost -args {z0,parameters} -lang:c++
codegen Wind_cost_GN -args {z0,parameters} -lang:c++

tic
for i=2:N_opt
    [zstar,fxstar,k,exitflag,~] = myfmincon(@(z)Wind_cost_mex(z,parameters),z0,A,b,C,d,p,q,myoptions);
    %z0                      = zstar;
    parameters.r_meas       = Position.signals.values(i,:)';
    parameters.rd_meas      = PositionDot.signals.values(i,:)';
    parameters.rdd_meas     = PositionDotDot.signals.values(i,:)';
    parameters.F_T_norm     = Forces.signals.values(i,end);
    parameters.Cl           = Cl.signals.values(i);
    parameters.Cd           = Cd.signals.values(i);
    W0_vec(i,:)             = zstar(1:3);
    parameters.zold         = zstar;
    myoptions.GN_funF       = @(z)Wind_cost_GN_mex(z,parameters);
    fprintf("Iteration %d done, Wind is [%f %f %f], Zl is [%f %f %f] (norm %f), iter: %d, exit:%d \n", ...
            i, zstar(1), zstar(2), zstar(3),zstar(4), zstar(5), zstar(6), norm([zstar(4) zstar(5) zstar(6)]), k, exitflag);
end
toc

figure(1)
plot(1:N_opt, W_log.signals.values(1:N_opt,1), 1:N_opt, W0_vec(:,1));
figure(2)
plot(1:N_opt, W_log.signals.values(1:N_opt,2), 1:N_opt, W0_vec(:,2));
%figure(3)
%plot(1:N_opt, W_log.signals.values(1:N_opt,3), 1:N_opt, W0_vec(:,3));
