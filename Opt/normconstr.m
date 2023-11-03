function [c,ceq,gc,gceq] = normconstr(z, deltaNorm)

    zsqrt = sqrt(z(4)^2+z(5)^2+z(6)^2);
    c = [zsqrt-1-deltaNorm;
         -zsqrt+1-deltaNorm]; 
    ceq = [];

    if nargout > 2 % gradient required
        gc = zeros(6,2);
        gc(4:6,1) = [z(4)/zsqrt;
                     z(5)/zsqrt;
                     z(6)/zsqrt];
        gc(4:6,2) = -gc(4:6,1);
        gceq = [];
    end
end
