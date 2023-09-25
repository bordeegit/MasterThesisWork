function u_rot = FLRot(u,v,psi)
    v = v / norm(v);
    
    % Calculate the skew-symmetric matrix of v
    skew_v = [0, -v(3), v(2);
              v(3), 0, -v(1);
              -v(2), v(1), 0];
    
    % Calculate the rotation matrix using the Rodrigues' rotation formula
    R = eye(3) + sin(psi) * skew_v + (1 - cos(psi)) * (skew_v^2);
    
    % Rotate vector u using the rotation matrix R
    u_rot = R * u;
end