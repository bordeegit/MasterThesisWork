function [xkp1] = StateEquations(xk, input, Ts, mk, mt_noL)
    
    r = xk(1:3);
    rdot = xk(4:6);
    rdotdot = xk(7:9);
    w_n = xk(10:11);
    F_L = xk(12:14);
    F_D_norm = xk(15);
    c_u = xk(16);

    u_s = input(1);
    F_T_norm = input(2);

    w_a = [w_n; 0] - rdot;   % w_a = w_n - rdot
    L = norm(r);
    g = [0;0;-9.81];
    m = mk + 0.25*mt_noL*L;

    xkp1 = zeros(16,1);

    xkp1(1:3) = r + rdot.*Ts;                          %r
    xkp1(4:6) = rdot + rdotdot.*Ts;                          %rd
    xkp1(7:9) = (F_L - F_T_norm*r/L + ...
                F_D_norm*w_a/norm(w_a) + (mk+0.5*mt_noL*L)*g)/m;  %rdd
    xkp1(10:11) = w_n;                                    %w_n
    xkp1(12:14) = FLRot(F_L,w_a,c_u*u_s*Ts);      %F_L vec
    xkp1(15) = F_D_norm;                                          %F_D
    xkp1(16) = c_u;                                          %cu
end

