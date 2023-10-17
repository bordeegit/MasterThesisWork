function [F] = Wind_cost(W,p)
   
    
    L = norm(p.r_meas);
    m = p.mk + 0.25*p.mt_noL*L;

    wa = W - p.rd_meas;
    wa_norm = norm(wa);
    Fg = (p.mk + 0.5*p.mt_noL*L)*[0;0;-p.g];
    Ft = p.F_T_norm*p.r_meas/L;
    % Q: Come gestire Fl? sto usando la misurazione ma non possiamo 
    % assumere venga misurato, però non conosciamo la direzione, magari con
    % lo steering input e facendo la rotazione?
    Fl = p.Fl;
    Fd = 0.5*p.rho*p.A*p.Cd*wa_norm^2*wa/wa_norm;
    
    rdd = (Fl + Fd + Fg - Ft)/m;

    F = norm(p.rdd_meas-rdd)^2;  % 2-norm squared

end

