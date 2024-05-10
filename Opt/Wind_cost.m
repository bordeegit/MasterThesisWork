function [ftot,gradtot] = Wind_cost(z,p)
   
    
    L = norm(p.r_meas);
    m = p.mk + 0.25*p.mt_noL*L;
    Fdconst = 0.5*p.rho*p.A*p.Cd;
    Flconst = 0.5*p.rho*p.A*p.Cl;

    wa = [z(1:2);0] - p.rd_meas;
    wa_norm = norm(wa);
    Fg = (p.mk + 0.5*p.mt_noL*L)*[0;0;-p.g];
    Ft = p.F_T_norm*p.r_meas/L;
    Fl = Flconst*wa_norm^2*z(3:5);
    Fd = Fdconst*wa_norm*wa;
    %Fd = Fdconst*wa_norm^2*wa/wa_norm;
    zdiff = z-p.zold;
    
    rdd = (Fl + Fd + Fg - Ft)/m;

    F = [sqrt(p.Q)*(p.rdd_meas-rdd);
         sqrt(p.Qw)*zdiff];
         %p.gamma*(z(3)*wa(1)+z(4)*wa(2)+z(5)*wa(3))];  % should be the same of 2-norm squared

    ftot = F'*F;

    if nargout > 1 % gradient required
        gradF = zeros(size(z,1),size(F,1));
        gradF(1,1) = -sqrt(p.Q(1,1))/m * (Flconst*z(3)*2*wa(1) + ...
            Fdconst*wa_norm + Fdconst*2*wa(1)^2/(2*wa_norm));
        gradF(2,1) = -sqrt(p.Q(1,1))/m * (Flconst*z(3)*2*wa(2) + ...
            Fdconst*2*wa(1)*wa(2)/(2*wa_norm));
        gradF(3,1) = -sqrt(p.Q(1,1))/m * Flconst*wa_norm^2;

        gradF(1,2) = -sqrt(p.Q(2,2))/m * (Flconst*z(4)*2*wa(1) + ...
            Fdconst*2*wa(2)*wa(1)/(2*wa_norm));
        gradF(2,2) = -sqrt(p.Q(2,2))/m * (Flconst*z(4)*2*wa(2) + ...
            Fdconst*wa_norm + Fdconst*2*wa(2)^2/(2*wa_norm));
        gradF(4,2) = -sqrt(p.Q(2,2))/m * Flconst*wa_norm^2;

        gradF(1,3) = -sqrt(p.Q(3,3))/m * (Flconst*z(5)*2*wa(1) + ...
            Fdconst*2*wa(3)*wa(1)/(2*wa_norm));
        gradF(2,3) = -sqrt(p.Q(3,3))/m * (Flconst*z(5)*2*wa(2) + ...
            Fdconst*2*wa(3)*wa(2)/(2*wa_norm));
        gradF(5,3) = -sqrt(p.Q(3,3))/m * Flconst*wa_norm^2;

        gradF(:,3:7) = sqrt(p.Qw);
        gradtot = 2*gradF*F;
    end
end

