function plotWindEstimation(windData, noZ)
    if ~exist('noZ', 'var')
        noZ = false;
    end
    
    n_subpl = 2 + (~noZ);
    components = {'x', 'y', 'z'};
    
    % Create figure for component-wise plots
    figure('Name', 'Wind Components');
    
    % Plot each component
    for i = 1:n_subpl
        subplot(n_subpl, 1, i);
        plotWindComponent(windData, i, components{i});
    end
    
    % Link x-axes
    linkaxes(findall(gcf, 'type', 'axes'), 'x');
    
    % Create magnitude and direction plots
    plotWindMagnitude(windData);
    plotWindDirection(windData);
end