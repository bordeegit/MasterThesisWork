function [y] = MeasurementEquations(x)
%MEASUREMENTEQUATIONS 
%   Function to compute the measurement vector y

    y = zeros(7,1);

    y(1:3) = x(1:3); % add noise maybe? TODO
    y(4:6) = x(4:6); % add noise maybe? 
    y(7) = x(12:14)'*([x(10:11); 0] - x(4:6));

end

