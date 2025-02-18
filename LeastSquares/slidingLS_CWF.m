function W_est_full = slidingLS_CWF(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ,CWF_mask)
%SLIDINGLS_CWF Solves the Estimation problem with a 'sliding window' approach only when a condition is met 
%   At each step, if the previous blockSize measurment are valid (considering CWF_mask)
%   computes the wind solving a least square problem considering the previous blockSize measurement

tol = 1e-15;
maxIter = 20;

iSequence = initStep+blockSize:maxStep;
A = zeros(blockSize,3-noZ);
B = zeros(blockSize,1);
W_est_vals = NaN(maxStep,3-noZ);
for i = iSequence
    if all(~CWF_mask(i-blockSize+1:i)) % Check if the block is valid 
        for j = 1:blockSize
            A(j,:) = l(i-blockSize+j,1:3-noZ); 
            B(j) = W_er_norm(i-blockSize+j)+L_dot(i-blockSize+j);
        end
        evalc("W_est_vals(i,:) = lsqr(A,B,tol,maxIter)';");
    end
end

if(noZ)
    W_est_vals = [W_est_vals zeros(length(W_est_vals),1)];
end

W_est_full = W_est_vals;

end