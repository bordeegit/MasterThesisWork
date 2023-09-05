function [xkp1] = StateEquations(xk, Ts)
%STATEEquations Summary of this function goes here
    
    F_L = 
    F_D = 
    m_k = mass;
    m_t = (1/4)*n_line*Line_diameter^2*Line_density*pi*norm(xk(1:3));
    A = inv([mass*eye(3), xk(1:3)
             xk(1:3)', 0]);
    B = 

    xkp1(1:3) = xk(1:3) + xk(4:6).*Ts;
    xkp1(4:6) = xk(4:6) + xk(7:9).*Ts;

end

