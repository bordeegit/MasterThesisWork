clear 
close all 

addpath(genpath(pwd));
load(which('Log_12m2_fix.mat'));
% Kite_phi/Kite_theta and phi/theta are basically the same, just 2
% different acquisitions, first from the ground, second from IMU (smoother)

% T_s = 0.02, acquisition lasted 500s -> 25000 steps + 1 @t=0

DataFlag = "RealFlight";
r = Kite_param.length;
N = length(theta);
T_s = 0.02;

parameters.rho = 1.2;
parameters.A = Kite_param.area;
parameters.Cd_l = 1.2;          %supposto da dati precedenti TODO
parameters.d_l = 0.003;         %supposto da dati precedenti TODO
parameters.n_l = Kite_param.nt;
parameters.mk = Kite_param.mass;
parameters.mt_noL = parameters.n_l*parameters.d_l^2*pi*970; 
parameters.g = 9.81;
parameters.Cd = 0.1516;
parameters.Cl = 0.85;

L_dot = zeros(N,1);



% Filter for Cable Force and Wind 
s=tf('s');
wb_filt=1;  % Hz
filt=1/(1+s/(wb_filt*2*pi));
filt=c2d(filt,T_s,'tustin');
[Bfilt,Afilt]=tfdata(filt,'v');

%% Kinematics

pos = r.*[cos(Kite_phi).*cos(Kite_theta),...
          sin(Kite_phi).*cos(Kite_theta),...
          sin(Kite_theta)];

% Velocity seems ok, the gradient version is a bit more noisy, it makes
% sense, but is in line with the analytical with the transformation R_GL

% To find this I can do the analytical derivative, or use R_gl*V_l, for
% which I already have the equations in the paper
posDot = zeros(N, 3);
for i=1:N
R_LG = [-cos(phi(i))*sin(theta(i)) -sin(phi(i)) -cos(phi(i))*cos(theta(i));
          -sin(phi(i))*sin(theta(i))  cos(phi(i)) -sin(phi(i))*cos(theta(i));
                cos(theta(i))          0           -sin(theta(i))   ];
posDot(i,:) = R_LG * (r*[theta_dot(i); cos(theta(i))*phi_dot(i); 0]);
end

% We don't need accelerationd, no need for posDotDot
% posDotDot = [0 0 0; diff(posDot)/T_s]; 


%% Wind

% Needs to be rotated, according to Log Analysis
W_unf = [Wind_speed.*cos(pi/2 - Wind_dir),...
         Wind_speed.*sin(pi/2 - Wind_dir).*(1+0.1*sin(theta)), ...
         zeros(N,1)];
W = filter(Bfilt,Afilt,W_unf);

%% Force

ScaleFactor = 500/3.5/0.22481;
Load_cell_left_filt=(filter(Bfilt,Afilt,Load_cell_left) - 0.5)*ScaleFactor;
Load_cell_right_filt=(filter(Bfilt,Afilt,Load_cell_right) - 0.5)*ScaleFactor;
Load_cell_center_filt=(filter(Bfilt,Afilt,Load_cell_center) - 0.5)*ScaleFactor*2;

F_T_norm = Load_cell_left_filt+Load_cell_right_filt+Load_cell_center_filt;


%% Aero

Cl_sim = parameters.Cl*ones(N,1);
Cd_sim = parameters.Cd*ones(N,1);

%% Polish

clearvars -except DataFlag r N T_s parameters L_dot pos posDot posDotDot...
    W F_T_norm Cl_sim Cd_sim phi theta phi_dot theta_dot


