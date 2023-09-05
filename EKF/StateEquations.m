function [xkp1] = StateEquations(xk, Ts, m)
%STATEEquations Summary of this function goes here
    xkp1 = zeros(15,1);
    A = [m*eye(3), xk(1:3)
             xk(1:3)', 0];
    B = [xk(13:15) 
        -xk(4:6)'*xk(1:3) + norm(xk(4:6))^2 + norm(xk(1:3))*norm(xk(4:6))];

    xkp1(1:3) = xk(1:3) + xk(4:6).*Ts;
    xkp1(4:6) = xk(4:6) + xk(7:9).*Ts;
    xkp1(7:10) = A\B;
    xkp1(11:12) = xk(11:12);
    xkp1(13:15) = xk(13:15);

end

