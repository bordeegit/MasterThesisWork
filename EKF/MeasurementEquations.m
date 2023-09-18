function [y] = MeasurementEquations(x, hr, h0)
%MEASUREMENTEQUATIONS 
%   Function to compute the measurement vector y
% Maybe r and rdot need to be computed with theta and phi, How can i pass
% them at each time step? the function cannot see the Workspace

    y = zeros(9,1);

    y(1:3) = x(1:3); % add noise maybe? TODO
    y(4:6) = x(4:6); % add noise maybe? 
    y(7) = log(hr/h0)/(log(x(3)/h0)) * norm(x(11:12));
    y(8) = atan(x(12)/x(11));
    y(9) = norm(x(1:3)*x(10));

end

