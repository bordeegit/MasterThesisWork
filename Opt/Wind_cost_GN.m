function [ftot] = Wind_cost_GN(z,p)
   
    
    L = norm(p.r_meas);
    m = p.mk + 0.25*p.mt_noL*L;

    wa = z(1:3) - p.rd_meas;
    wa_norm = norm(wa);
    Fg = (p.mk + 0.5*p.mt_noL*L)*[0;0;-p.g];
    Ft = p.F_T_norm*p.r_meas/L;
    Fl = 0.5*p.rho*p.A*p.Cl*wa_norm^2*z(4:6);
    Fd = 0.5*p.rho*p.A*p.Cd*wa_norm^2*wa/wa_norm;
    
    rdd = (Fl + Fd + Fg - Ft)/m;

    F = [p.Q*(p.rdd_meas-rdd);
         p.Qw*(z-p.zold)];  % should be the same of 2-norm squared

    %g = [z(4)*(z(1)-p.rd_meas(1))+z(5)*(z(2)-p.rd_meas(2))+z(6)*(z(3)-p.rd_meas(3))];
         %sqrt(z(4)^2+z(5)^2+z(6)^2)-1];
    
    h = [-(z(4)*wa(1)+z(5)*wa(2)+z(6)*wa(3))+p.deltaOrth;
         (z(4)*wa(1)+z(5)*wa(2)+z(6)*wa(3))+p.deltaOrth;
         -sqrt(z(4)^2+z(5)^2+z(6)^2)+1+p.deltaNorm;
         sqrt(z(4)^2+z(5)^2+z(6)^2)-1+p.deltaNorm];
    %g = 0;

    ftot = [F;h];
end

