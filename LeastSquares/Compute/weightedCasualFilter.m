function filtered_wind = weightedCasualFilter(wind_data, window_size, prev_weight)
    % wind_data: input data array
    % window_size: size of the moving window
    % prev_weight: weight given to previous filtered value (0-1)
    
    filtered_wind = zeros(size(wind_data));
    filtered_wind(1) = wind_data(1);
    
    for i = 2:length(wind_data)
        % Get window of previous raw values
        if i < window_size
            window = wind_data(1:i);
        else
            window = wind_data(i-window_size+1:i);
        end
        
        % Remove NaN from current window
        valid_values = window(~isnan(window));
        
        if ~isempty(valid_values)
            % Current median from valid values
            current_median = median(valid_values);
            
            % Check if coming from NaN region
            if isnan(wind_data(i-1))
                % Weighted combination with previous filtered value
                filtered_wind(i) = prev_weight * filtered_wind(i-1) + ...
                                 (1-prev_weight) * current_median;
            else
                % Normal case - still use some weighting but less
                filtered_wind(i) = (prev_weight/2) * filtered_wind(i-1) + ...
                                 (1-prev_weight/2) * current_median;
            end
        else
            % If no valid values, maintain last filtered value
            filtered_wind(i) = filtered_wind(i-1);
        end
    end
end