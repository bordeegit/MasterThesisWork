function [xkp1] = StateEquations(xk, input, Ts, mk, mt_noL)
    
    w_a = [xk(11:12); 0] - xk(4:6);   % w_a = w_n - rdot
    L = norm(xk(1:3));
    g = [0;0;-9.81];
    m = mk + 0.25*mt_noL*L;

    xkp1 = zeros(16,1);

    xkp1(1:3) = xk(1:3) + xk(4:6).*Ts;                          %r
    xkp1(4:6) = xk(4:6) + xk(7:9).*Ts;                          %rd
    xkp1(7:9) = (xk(12:14) - input(2)*xk(1:3)/L + ...
                xk(16)*w_a/norm(w_a) + (mk+0.5*mt_noL*L)*g)/m;                                           %rdd, v
    xkp1(10:11) = xk(10:11);                                    %w_n
    xkp1(12:14) = FLRot(xk(12:14),w_a,xk(16)*input(1)*Ts);          %F_L vec
    xkp1(15) = xk(15);                                          %F_D
    xkp1(16) = xk(16);                                          %cu
end

