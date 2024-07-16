function [W_est_full, iSequence] = blockLS(W_er_norm,l,L_dot,blockSize,initStep,maxStep,noZ)
%BLOCKLS Solves the Estimation problem with a 'block' approach
%   Computes the wind speed at i considering the previous i-blockSize 
%   measurments; the next estimation will take place at 
%   i+blockSize, considering again the previous blockSize block
%   Outputs the full vector W_est_full, the values computed every blockSize
%   steps and the sequence of intervals at which the values where computed

tol = 1e-15;
maxIter = 20;

iSequence = initStep+blockSize:blockSize:maxStep;%-mod(maxStep-initStep,blockSize);
A = zeros(blockSize,3-noZ);
B = zeros(blockSize,1);
W_est_vals = zeros(length(iSequence),3-noZ);
for idx = 1:length(iSequence)
    i = iSequence(idx);
    for j = 1:blockSize
        A(j,:) = l(i-blockSize+j,1:3-noZ); 
        B(j) = W_er_norm(i-blockSize+j)+L_dot(i-blockSize+j);
    end
    W_est_vals(idx,:) = lsqr(A,B,tol,maxIter)';
    % Tried to use lsqminnorm, the same as pinv(A)*B, get same results to 10e-8
end

if(noZ)
    W_est_vals = [W_est_vals zeros(length(W_est_vals),1)];
end

W_est_full = [NaN(iSequence(1)-1,3); kron(W_est_vals, ones(blockSize,1))]; % Reconstructing full vector

end