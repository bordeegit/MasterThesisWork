function [c,ceq,gc,gceq] = normconstr(z, rd_meas)

    wa = [z(1:2);0] - rd_meas;
    zsqrt = sqrt(z(3)^2+z(4)^2+z(5)^2);
    zorth = z(3)*wa(1)+z(4)*wa(2)+z(5)*wa(3);
    c = []; 
    ceq = [zsqrt-1;
           zorth];

    if nargout > 2 % gradient required
        gc = [];
        gceq = zeros(5,2);
        gceq(3:5,1) = [z(3)/zsqrt;
                       z(4)/zsqrt;
                       z(5)/zsqrt];
        gceq(:,2) = [z(3:4);
                     wa];
    end
end
