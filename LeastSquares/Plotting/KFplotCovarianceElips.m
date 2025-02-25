function KFplotCovarianceElips(P)

N = length(P);

% Assuming P is your 2x2xN covariance matrix
figure;
hold on;

% Choose how many steps to visualize (can get cluttered if all N)
steps_to_show = min(N, 20); % Show up to 20 steps
step_size = max(1, floor(N/steps_to_show));

% Create colormap from light to dark
colors = cool(steps_to_show);

% Plot ellipses at selected timesteps
for i = 1:step_size:N
    step_idx = ceil(i/step_size);
    if step_idx <= steps_to_show
        % Get current 2x2 covariance
        current_P = squeeze(P(:,:,i));
        
        % Calculate ellipse points (95% confidence)
        [eigvec, eigval] = eig(current_P);
        chi_val = 5.991; % 95% confidence for 2 DOF
        theta = linspace(0, 2*pi, 100);
        
        % Generate ellipse points
        a = sqrt(chi_val * eigval(1,1));
        b = sqrt(chi_val * eigval(2,2));
        ellipse = [a*cos(theta); b*sin(theta)];
        ellipse = eigvec * ellipse;
        
        % Plot with color indicating timestep
        plot(ellipse(1,:), ellipse(2,:), 'Color', colors(step_idx,:), 'LineWidth', 1.5);
        
        % Add timestep label
        text(ellipse(1,1), ellipse(2,1), num2str(i), 'Color', colors(step_idx,:));
    end
end

title('Evolution of Covariance Ellipses');
xlabel('State Dimension 1');
ylabel('State Dimension 2');
grid on;
axis equal;
colormap(cool);
colorbar('Ticks', linspace(0,1,5), 'TickLabels', {'Start', '', 'Time', '', 'End'});