function plotWindDirection(windData)
    figure('Name', 'Wind Direction');
    grid on;
    hold on;
    
    fields = {'estFilt', 'CWF', 'OPT', 'real'};
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