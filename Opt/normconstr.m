function [c,ceq,gc,gceq] = normconstr(z, rd_meas)

    wa = [z(1:2);0] - rd_meas;
    zsqrt = sqrt(z(3)^2+z(4)^2+(1-z(3)^2-z(4)^2)^2); % simplify?
    zorth = z(3)*wa(1)+z(4)*wa(2)+(1-z(3)^2-z(4)^2)*wa(3);
    c = []; 
    ceq = [zsqrt-1;
           zorth];

    if nargout > 2 % gradient required
        gc = [];
        gceq = zeros(4,2);
        gceq(3:4,1) = [(4*z(3)^3+4*z(3)*z(4)^2-2*z(3))/(2*zsqrt);
                       (4*z(4)^3+4*z(4)*z(3)^2-2*z(4))/(2*zsqrt)];
        gceq(:,2) = [z(3:4);
                     wa(1)-2*z(3)*wa(3);
                     wa(2)-2*z(4)*wa(3)];
    end
end
