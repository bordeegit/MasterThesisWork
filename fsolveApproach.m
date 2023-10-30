close all 


%ind = 427 first time we arrive at F_T > 10000
ind = 500;
W_meas               = W_log.signals.values(ind,:)';
pm.r_meas            = Position.signals.values(ind,:)';
pm.rd_meas           = PositionDot.signals.values(ind,:)';
pm.F_T_norm_meas     = Forces.signals.values(ind,end);
pm.Cl_meas           = Cl.signals.values(ind);
pm.Cd_meas           = Cd.signals.values(ind);
pm.rdd_meas          = PositionDotDot.signals.values(ind,:)';
pm.mk                = mass;
pm.A                 = area;
pm.rho               = rho;
pm.mt_noL            = mt_noL;
pm.Fl_meas           = Forces.signals.values(ind,4:6)';
pm.zl_meas           = Fl_meas/norm(Fl_meas); 

 

fun = @(x)accDiff(x,pm);
%x0 = W_meas; 
x0 = [0;0;0];
x = fsolve(fun, x0)

fprintf("Estimated Wind is [%.2f,%.2f,%.2f]\nActual Wind is [%.2f,%.2f,%.2f]\n", ...
    x(1), x(2), x(3), W_meas(1),W_meas(2),W_meas(3));

function F = accDiff(x, pm)

    x(3) = 0;
    L = norm(pm.r_meas);
    m = pm.mk + 0.25*pm.mt_noL*L;
    wa = [x(1);x(2);x(3)] - pm.rd_meas;
    wa_norm = norm(wa);
    Fg = (pm.mk + 0.5*pm.mt_noL*L)*[0;0;-9.81];
    Ft = pm.F_T_norm_meas*pm.r_meas/L;
    Fl = 0.5*pm.rho*pm.A*pm.Cl_meas*wa_norm^2*pm.zl_meas;
    Fd = 0.5*pm.rho*pm.A*pm.Cd_meas*wa_norm^2*wa/wa_norm;

    rdd = (Fl + Fd + Fg - Ft)/m;

    F = rdd - pm.rdd_meas; 
end
