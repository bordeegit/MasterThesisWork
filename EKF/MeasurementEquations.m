function [y] = MeasurementEquations(x)

    y = zeros(6,1);

    y(1:3) = x(1:3); 
    y(4:6) = x(4:6); 

end

