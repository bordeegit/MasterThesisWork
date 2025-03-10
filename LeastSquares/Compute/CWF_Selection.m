function [W_est_CWF, CWF_mask, plotData] = CWF_Selection(W_est, phi_dot, CWF_TH_UPP, CWF_TH_LOW)

maxStep = length(W_est);
CWF_mask = phi_dot(2:maxStep+1)<= CWF_TH_UPP & phi_dot(2:maxStep+1)>= CWF_TH_LOW;
phi_dot_masked = phi_dot(2:maxStep);
phi_dot_masked(CWF_mask) = NaN;
W_est_CWF = W_est;
W_est_CWF(CWF_mask,:) = NaN;

plotData.original_phi_dot = phi_dot(2:maxStep);
plotData.masked_phi_dot = phi_dot_masked;
plotData.plot_options_original = {'Color', [0.565 0.808 0.98 0.2]};
plotData.plot_options_masked = {'Color', [0 0.4470 0.7410], 'LineWidth', 2};
plotData.title = '$\dot{\phi}$ masking';

if nargout < 3
    figure, grid on, box on, 
    plot( phi_dot(2:maxStep), 'Color', [0.565 0.808 0.98 0.2]), hold on
    plot(phi_dot_masked, 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
    title('$\dot{\phi}$ masking')
end
