% Init. file for model Kite_sim.mdl
% L. Fagiano
% fagiano@control.ee.ethz.ch

clear 
%close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Sampling time %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts_10ms=1e-2;                                   % 10 ms sampling time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Controller's param. %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_fb_K=2*2.8*0.01;                          % m/rad error
HLC_input_max=0.25;                             % m
Targets_matrix=[0.3 -0.3;
                0.3 0.3];                       % theta phi
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Kite & simulation param. %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim_stop_time=120;                              % s simulation stop time

rho=1.2;                                        % Air density (kg/m^3)
g=9.81;                                         % Gravity (m/s^2)
                                                
                                                % theta = 0 is vertical
                                                % theta = 90 is horizontal
                                                % As paper
init_stat  = [60*pi/180 1*pi/180 ...            % theta, phi
    1*pi/180 1*pi/180 ...                       % thetadot, phidot
    50 0];                                      % r, rdot

HLC2dl=-4;                                      % m actuator to m line diff.

area=9;                                         % Kite area (m^2)
mass=2.4;                                       % Kite mass (kg)
wingspan=2;                                     % Kite wing span (m)
    
alpha_0=2.5*pi/180;                             % Base AoA
alpha_min=-8*pi/180;
alpha_max=18*pi/180;

CD_Line=1.2;                                    % Coeff. drag tether
Line_diameter=0.003;                            % Tether diam.
Line_density=970;                               % Tether density (kg/m^3)
n_line=3;                                       % n. of lines
mt_noL = n_line*Line_diameter^2*pi*Line_density; 

%%%%% Wind Parameters %%%%%
%% SHEAR WIND MODEL
%     hr=32.5;                                        % Wind shear model
%     phi_wr = 30*pi/180;
%     %phi_wr = 150*pi/180;
%     h0=6e-4; %BASE                                 % surface roughness
%     %h0 = 2;
%     %In Sheat Wind Estimation paper we get h0 =
%     %   0.15 for Category C flight phases(terminal flight phases, which 
%     %                                     include takeoff, approach,and
%     %                                     landing)
%     %   2.0 otherwise
%     wr=4;

%% WIND MODEL

%%% Constant Wind %%%
ConstWindX = 5;                                      % Wind dist. X
ConstWindY = 0;                                      % Wind dist. Y
ConstWindZ = 0;

%%% WRT XY %%%
x_pos_Data = [0,50,30,60,60,60,25,0,0,100,100,100,80];
y_pos_Data = [0,10,-10,0,-50,50,0,-50,50,0,-50,50,-20];
w_y_pos_Data = [0,10,5,0,0,0,0,0,0,0,0,0,-10];


%%% WRT Z %%%
z_scale = 0.7;
height_Data = [0 5 10 15 20 25 40 50];
w_x_height_Data = [0 2 8 10 12 20 20 20];

%% Curve Fitting Approach (Slow and useless)
% [xDataCurve, yDataCurve] = prepareCurveData( height_Data, w_x_height_Data );
% 
% % Set up fittype and options
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );
% opts.SmoothingParam = 0.8;
% 
% % Fit model to data.
% [WxFunctionZ, ~] = fit( xDataCurve, yDataCurve, ft, opts );
% 
% %%% WRT XY %%%
% x_pos_Data = [0,50,30,60,60,60,25,0,0,100,100,100,80];
% y_pos_Data = [0,10,-10,0,-50,50,0,-50,50,0,-50,50,-20];
% w_y_pos_Data = [0,10,5,0,0,0,0,0,0,0,0,0,-10];
% 
% [xDataSurf, yDataSurf, zDataSurf] = prepareSurfaceData( x_pos_Data, y_pos_Data, w_y_pos_Data );
% 
% % Fit model to data.
% [WyFunctionXY, ~] = fit( [xDataSurf, yDataSurf], zDataSurf, 'cubicinterp', 'Normalize', 'on' );
% 
% % Creating Function handles for simulink
% WxFncH = @(z)feval(WxFunctionZ,z);
% WyFncH = @(xy)feval(WyFunctionXY,xy);


%% CONTROL
K_r= 1000;                                      % reeling control gain
r_ref=init_stat(5);                             % reeling ref.
r_slope=0;                                      % reeling ref.

alpha_var=[-8.5 0 6 10 14 16 18 20]*pi/180';    % AoA table (deg)
CL_var=[-0.147 0.332 0.652 0.838...
    0.99 1.052 1.084 1.032]'*0.8;               % Lift coeff. table
CD_var=[0.010 0.013 0.048 0.083...
    0.128 0.154 0.181 0.218]';                  % Drag coeff. table


