function [W_est_full, iSequence] = slidingLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ)
%SLIDINGLS Solves the Estimation problem with a 'sliding window' approach 
%   At each step, computes the wind solving a least square problem
%   considering the previous blockSize measurement

tol = 1e-15;
maxIter = 20;

iSequence = initStep+blockSize:maxStep;
A = zeros(blockSize,3-noZ);
B = zeros(blockSize,1);
W_est_vals = [NaN(iSequence(1)-1,3-noZ); zeros(length(iSequence),3-noZ)];
for i = iSequence
    for j = 1:blockSize
        A(j,:) = l(i-blockSize+j,1:3-noZ); 
        B(j) = W_er_norm(i-blockSize+j)+L_dot(i-blockSize+j);
    end
    W_est_vals(i,:) = lsqr(A,B,tol,maxIter)';
    % Tried to use lsqminnorm, the same as pinv(A)*B, get same results to 10e-8
end

if(noZ)
    W_est_vals = [W_est_vals zeros(length(W_est_vals),1)];
end

W_est_full = W_est_vals;

end