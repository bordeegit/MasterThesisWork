function [x_dot]=kiteEq(x_t,Psi,F_line,W,p)

% Discrete-time kite dynamics with fixed KSU (e.g. for
% KE-yoyo configuration)

% Unpacking parameters

A           = p.A;
m           = p.m;
rho         = p.rho;
g           = p.g;
n_l         = p.n_l;
d_l         = p.d_l;
CD_l        = p.CD_l;
rho_l       = p.rho_l;
alpha_0     = p.alpha_0;
alpha_max   = p.alpha_max;
alpha_min   = p.alpha_min;
alpha_var   = p.alpha_var;
Cl_var      = p.Cl_var;
Cd_var      = p.Cd_var;

% Unpacking states

theta       =       x_t(1);
phi         =       x_t(2);
r           =       x_t(5);
thetadot    =       x_t(3);     
phidot      =       x_t(4);
rdot        =       x_t(6);
psi         =       Psi;
                                                                                                                      
ep=[...
    cos(theta)*cos(phi)  -sin(phi)    sin(theta)*cos(phi)
    cos(theta)*sin(phi)  cos(phi)     sin(theta)*sin(phi)
    -sin(theta)          0            cos(theta)
    ];                               % local coordinate system -> nominal wind coordinate system


wp=ep'*W;                           % wind in local coordinate system

vp=[...
    r*thetadot
    r*sin(theta)*phidot
    rdot
    ];                                % apparent wind due to kite movement (local system)
vp=-vp;

wep=vp+wp;                          % calculate effective wind we (in local coordinates)

weo=ep*wep;                         % effective wind we in nominal wind coordinate system
weobar=norm(weo);                   % norm of effective wind speed

epp=zeros(3,3);

epp(:,1)=wep;                       % epp(:,1) pointing in projected direction of effective wind wep
epp(3,1)=0;
wbar=norm(epp(:,1));
wbar=max(1e-8,wbar);                % saturation on minimal wind tangential to the sphere

epp(:,1)=epp(:,1)/wbar;             % projected coordinate system - axis 1

epp(:,3)=[0;0;1];                   % projected coordinate system - axis 3

epp(:,2)=[-epp(2,1),epp(1,1),0];    % projected coordinate system - axis 2

wepp=epp'*wep;                      % translate wep in new system

ratio=wepp(3)*tan(psi)/wepp(1);

ratio=min(1,ratio);                 % saturation on sine of eta angle
ratio=max(-1,ratio);                % saturation on sine of eta angle

eta = asin(ratio);

etpp=[cos(psi)*(-sin(eta));cos(psi)*(cos(eta));sin(psi)];

edpp=wepp/norm(wepp);               % drag direction
ldpp=outprod(edpp,etpp);            % lift direction

edp=epp*edpp;                       % Turn edpp into local system
ldp=epp*ldpp;                       % Turn ldpp into local system

Pdyn=weobar*weobar*1/2*rho;

tan_d_alpha=wepp(3,1)/wepp(1,1);
d_alpha=atan(tan_d_alpha);          % angle of attack variation

alpha=alpha_0+d_alpha;              % angle of attack
alpha=min(alpha_max,alpha);         % angle of attack saturation
alpha=max(alpha_min,alpha);         % angle of attack saturation

Cl=interp1(alpha_var,Cl_var,alpha,'linear','extrap'); % lift coefficient
Cd=interp1(alpha_var,Cd_var,alpha,'linear','extrap'); % drag coefficient


Fl=Pdyn*A*Cl*ldp;                   % lift force

Drag_line=Pdyn*1/4*cos(d_alpha)*r*d_l*n_l*CD_l;     % Cable drag force
Fd_tot=(Pdyn*A*Cd+Drag_line)*edp;                   % Total drag force

Fapp=m*[phidot^2*r*sin(theta)*cos(theta)-2*rdot*thetadot;
    -2*rdot*phidot*sin(theta) - 2*thetadot*phidot*r*cos(theta);
    r*thetadot^2+r*phidot^2*sin(theta)^2]; % apparent forces

Fgrav=(m+n_l*rho_l*pi*r*d_l^2/4)*g*[sin(theta);0;-cos(theta)]; % gravity forces

% F_cavo=[0;0;(Fl(3,1)+Fapp(3,1)+Fd_tot(3,1)+Fgrav(3,1))-K_r*(r_slope-rdot)];
% F_cavo(3,1)=max(0,F_cavo(3,1)); % cable traction force (always positive and subtracted to other forces)

F_tot=Fl+Fd_tot+Fapp+Fgrav-F_line;

thetadotdot = F_tot(1,1)/(m*r);
phidotdot = F_tot(2,1)/(m*r*sin(theta));
rdotdot = F_tot(3,1)/m;

x_dot=[ thetadot;                           % theta dot
        phidot;                             % phi dot
        thetadotdot;                        % theta dotdot
        phidotdot;                          % phi dotdot
        rdot;                               % r dot
        rdotdot];                           % r dotdot

function c=outprod(a, b)
c=double(zeros(3,1));
c(1)=  a(2)*b(3) - a(3)*b(2);
c(2)=  a(3)*b(1) - a(1)*b(3);
c(3)=  a(1)*b(2) - a(2)*b(1);
c=c(:);


