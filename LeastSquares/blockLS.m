function [W_est_rec, W_est_vals, iSequence] = blockLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ)
    
    tol = 1e-15;
    maxIter = 20;
    
    iSequence = initStep:blockSize:maxStep-initStep;
    A = zeros(blockSize,3-noZ);
    B = zeros(blockSize,1);
    W_est_vals = zeros(floor((maxStep)/blockSize),3-noZ);
    for idx = 1:length(iSequence)
        i = iSequence(idx);
        for j = 1:blockSize
            A(j,:) = l(i+j-1,1:3-noZ); 
            B(j) = W_er_norm(i+j-1)+L_dot(i+j-1);
        end
        W_est_vals(idx,:) = lsqr(A,B,tol,maxIter)';
        % Tried to use lsqminnorm, the same as pinv(A)*B, get same results to 10e-8
    end
    
    if(noZ)
        W_est_vals = [W_est_vals zeros(length(W_est_vals),1)];
    end

    W_est_rec = kron(W_est_vals, ones(blockSize,1)); % Reconstructing full vector

end