function [x,P]=StepEKF(fstate,x,P,hmeas,y,Q,R, diffTs, us)

    [x1,F] = GradientF(fstate,x, us,'FD',diffTs);    %nonlinear update and linearization at current state
    F = F';
    P=F*(P*F) + Q;                 %partial update
    [y1,H] = GradientH(hmeas,x1,'FD',diffTs);    %nonlinear measurement and linearization
    H = H';
    zhat = y - y1;
    S = H*P*H' + R;
    K = P*H'*inv(S);
    x = x1 + K*zhat;
    P = (eye(size(x,1))-K*H)*P;
    
    
%     P12=P*H';                   %cross covariance
%     % K=P12*inv(H*P12+R);       %Kalman filter gain
%     % x=x1+K*(z-z1);            %state estimate
%     % P=P-K*P12';               %state covariance matrix
%     R=chol(H*P12+R);            %Cholesky factorization
%     U=P12/R;                    %K=U/R'; Faster because of back substitution
%     x=x1+U*(R'\(z-z1));         %Back substitution to get state update
%     P=P-U*U';                   %Covariance update, U*U'=P12/R/R'*P12'=K*P12.

end

