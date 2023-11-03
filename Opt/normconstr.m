function [c,ceq] = normconstr(z, deltaNorm)

    c = [sqrt(z(4)^2+z(5)^2+z(6)^2)-1-deltaNorm;
         -sqrt(z(4)^2+z(5)^2+z(6)^2)+1-deltaNorm]; 
    ceq = [];

end
