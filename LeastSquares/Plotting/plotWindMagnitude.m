function plotWindMagnitude(windData, tabgroup)
    if isnumeric(tabgroup)
        figure('Name', 'Wind Magnitude');
    else 
        tab = uitab(tabgroup, 'Title', 'Wind Magnitude');
        axes('Parent', tab);
    end
    
    grid on;
    hold on;
    
    fields = windData.plotFields;
    legendEntries = {};
    
    for i = 1:length(fields)
        field = fields{i};
        data = windData.estimates.(field);
        props = windData.plotProps.(field);
        
        magnitude = vecnorm(data')';
        plot(windData.time, magnitude, props);
        legendEntries{end+1} = sprintf('$|W_{%s}|$', field);
    end
    
    legend(legendEntries, 'Interpreter', 'latex');
    ylabel('Wind Speed (m/s)');
    xlabel('Time (s)');
    hold off;
end