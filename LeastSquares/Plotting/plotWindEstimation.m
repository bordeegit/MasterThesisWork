function plotWindEstimation(windData, noZ, tabgroup)
    if ~exist('noZ', 'var')
        noZ = false;
    end
    
    n_subpl = 2 + (~noZ);
    components = {'x', 'y', 'z'};
    
    % Create figure for component-wise plots
    if nargin < 3
        figure('Name', 'Wind Components');
        tabgroup = NaN;
    else 
        tab = uitab(tabgroup, 'Title', 'Wind Estimation Details');
        axes('Parent', tab);
    end
    
    sbpHandles = gobjects(n_subpl,1);
    % Plot each component
    for i = 1:n_subpl
        sbpHandles(i) = subplot(n_subpl, 1, i);
        plotWindComponent(windData, i, components{i});
    end
    
    % Link x-axes
    linkaxes(sbpHandles, 'x');
    
    % Create magnitude and direction plots
    plotWindMagnitude(windData, tabgroup);
    plotWindDirection(windData, tabgroup);
end