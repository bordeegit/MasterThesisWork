function plotWindComponent(windData, component, label)
    grid on;
    hold on;
    
    fields = windData.plotFields;
    legendEntries = {};
    
    for i = 1:length(fields)
        field = fields{i};
        data = windData.estimates.(field);
        props = windData.plotProps.(field);
        
        plot(windData.time, data(:,component), props);
        legendEntries{end+1} = sprintf('$W_%s$ (%s)', label, field);
    end
    
    legend(legendEntries, 'Interpreter', 'latex');
    %ylim([-5 7]);
    ylabel('Wind Speed (m/s)');
    xlabel('Time (s)');
    hold off;
end