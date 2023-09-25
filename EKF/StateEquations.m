function [xkp1] = StateEquations(xk, u_sk, Ts, mk, mt_noL)
%STATEEquations Summary of this function goes here
    
    w_a = [xk(11:12); 0] - xk(4:6);   % w_a = w_n - rdot
    L = norm(xk(1:3));
    g = [0,0,9.81]';
    m = mk + mt_noL*L; % this may not be correct, ask

    xkp1 = zeros(17,1);
    A = [m*eye(3), xk(1:3) ; xk(1:3)', 0];
    B = [xk(13:15) + xk(17)*w_a/norm(w_a) + (mk+0.5*mt_noL*L)*g; 
        -xk(4:6)'*xk(1:3) + norm(xk(4:6))^2 + L*norm(xk(7:9))];

    xkp1(1:3) = xk(1:3) + xk(4:6).*Ts;                          %r
    xkp1(4:6) = xk(4:6) + xk(7:9).*Ts;                          %rd
    xkp1(7:10) = A\B;                                           %rdd, v
    xkp1(11:12) = xk(11:12);                                    %w_n
    xkp1(13:15) = FLRot(xk(13:15),w_a,xk(17)*u_sk*Ts);          %F_L vec
    xkp1(16) = xk(16);                                          %F_D
    xkp1(17) = xk(17);                                          %cu
end

