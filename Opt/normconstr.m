function [c,ceq,gc,gceq] = normconstr(z, rd_meas)

    wa = z(1:3) - rd_meas;
    zsqrt = sqrt(z(4)^2+z(5)^2+z(6)^2);
    zorth = z(4)*wa(1)+z(5)*wa(2)+z(6)*wa(3);
    c = []; 
    ceq = [zsqrt-1;
           zorth];

    if nargout > 2 % gradient required
        gc = [];
        gceq = zeros(6,2);
        gceq(4:6,1) = [z(4)/zsqrt;
                       z(5)/zsqrt;
                       z(6)/zsqrt];
        gceq(:,2) = [z(4:6);
                     wa];
    end
end
