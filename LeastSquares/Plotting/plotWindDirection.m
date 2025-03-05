function plotWindDirection(windData, tabgroup)
    if isnumeric(tabgroup)
        figure('Name', 'Wind Direction');
    else 
        tab = uitab(tabgroup, 'Title', 'Wind Direction');
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
        
        direction = atan2d(data(:,2), data(:,1));
        plot(windData.time, direction, props);
        legendEntries{end+1} = sprintf('$\\angle W_{%s}$', field);
    end
    
    legend(legendEntries, 'Interpreter', 'latex');
    ylabel('Wind Direction (deg)');
    xlabel('Time (s)');
    hold off;
end