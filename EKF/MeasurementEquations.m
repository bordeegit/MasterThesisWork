function [y] = MeasurementEquations(x)

    y = zeros(7,1);

    y(1:3) = x(1:3); 
    y(4:6) = x(4:6); 
    y(7) = x(12:14)'*([x(10:11); 0] - x(4:6));

end

