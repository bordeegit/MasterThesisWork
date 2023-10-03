function [y] = MeasurementEquations(x)
%MEASUREMENTEQUATIONS 
%   Function to compute the measurement vector y

    y = zeros(8,1);

    y(1:3) = x(1:3); % add noise maybe? TODO
    y(4:6) = x(4:6); % add noise maybe? 
    y(7) = norm(x(1:3)*x(10));
    y(8) = x(13:15)'*([x(11:12); 0] - x(4:6));

end

