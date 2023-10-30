function [ftot] = Wind_cost(z,p)
   
    
    L = norm(p.r_meas);
    m = p.mk + 0.25*p.mt_noL*L;

    wa = z(1:3) - p.rd_meas;
    wa_norm = norm(wa);
    Fg = (p.mk + 0.5*p.mt_noL*L)*[0;0;-p.g];
    Ft = p.F_T_norm*p.r_meas/L;
    % Q: Come gestire Fl? sto usando la misurazione ma non possiamo 
    % assumere venga misurato, però non conosciamo la direzione, magari con
    % lo steering input e facendo la rotazione?
    Fl = 0.5*p.rho*p.A*p.Cl*wa_norm^2*z(4:6);
    Fd = 0.5*p.rho*p.A*p.Cd*wa_norm^2*wa/wa_norm;
    
    rdd = (Fl + Fd + Fg - Ft)/m;

    f = [norm(p.rdd_meas-rdd)^2
        p.Qw*(z-p.zold)];  % 2-norm squared

    g = z(4)*(z(1) - p.rd_meas(1)) + z(5)*(z(2) - p.rd_meas(2)) + z(6)*(z(3) - p.rd_meas(3));
    
    ftot = [f;g];
end

