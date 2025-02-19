function PerformanceEvaluation(W_matrices, PerfNames)

W_real = W_matrices{end};
PerfRMSE = zeros(length(PerfNames), 4);
PerfMean = zeros(length(PerfNames), 3);
ind = 1; 
for i = 1:length(W_matrices)
    W = W_matrices{i}; 
    PerfRMSE(ind,1:3) = rmse(W_real, W, 'omitnan');
    PerfRMSE(ind,4) = norm(PerfRMSE(ind,1:3));
    PerfMean(ind,:) = mean(W, 1, 'omitnan'); % Mean along each column
    ind = ind + 1;
end
PerfTable = table( PerfRMSE(:,1), PerfRMSE(:,2), PerfRMSE(:,4), PerfMean(:,1), PerfMean(:,2), ...
    'VariableNames', [ "RMSEx", "RMSEy", "RMSE_norm", "MeanX", "MeanY"], ...
    'RowNames',PerfNames);
disp(PerfTable)

end 