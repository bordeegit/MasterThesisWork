function [ftot] = Wind_cost(z,p)
   
    
    L = norm(p.r_meas);
    m = p.mk + 0.25*p.mt_noL*L;

    wa = z(1:3) - p.rd_meas;
    wa_norm = norm(wa);
    Fg = (p.mk + 0.5*p.mt_noL*L)*[0;0;-p.g];
    Ft = p.F_T_norm*p.r_meas/L;
    Fl = 0.5*p.rho*p.A*p.Cl*wa_norm^2*z(4:6);
    Fd = 0.5*p.rho*p.A*p.Cd*wa_norm^2*wa/wa_norm;
    zdiff = z-p.zold;
    
    rdd = (Fl + Fd + Fg - Ft)/m;

    F = [p.Q*(p.rdd_meas-rdd);
         p.Qw*zdiff];
         %p.gamma*(z(4)*wa(1)+z(5)*wa(2)+z(6)*wa(3))];  % should be the same of 2-norm squared

    ftot = F'*F;
end

