function printTrajNaN(pos, posDot, N_start, N_end, color, mask)
    %figure;
    grid on; hold on;
    
    % Determine colormap and colorbar label
    if nargin > 4 && ~isempty(color)
        colormap(color);
        barlabel = 'Unfilt. Est. Error';
    else
        colormap(winter);
        barlabel = 'Ground Speed (m/s)';
    end

    speed = vecnorm(posDot, 2, 2); % Compute speed magnitude

    % Use mask if provided; otherwise, use default range
    if nargin == 6 && ~isempty(mask)
        validIndices = find(~mask); % Get indices where mask is true (valid points)
    else
        validIndices = N_start:N_end;
    end

    % Break trajectory into segments to avoid drawing through masked points
    nanIdx = [diff(validIndices) > 1; true]; % Detect gaps in indices
    segmentStart = [1; find(nanIdx) + 1];
    segmentEnd = [find(nanIdx); numel(validIndices)];

    % Plot each valid segment separately
    for i = 1:numel(segmentStart)
        idxRange = validIndices(segmentStart(i):segmentEnd(i));
        if ~isempty(idxRange)
            patch([pos(idxRange,1); nan], ...
                  [pos(idxRange,2); nan], ...
                  [pos(idxRange,3); nan], ...
                  [speed(idxRange); nan], ...
                  'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.5);
        end
    end

    % Colorbar settings
    cb = colorbar;
    ylabel(cb, barlabel, 'Rotation', 270);

    % 3D view settings
    view(3);
    plot3(0, 0, 0, 'k*'); % Origin marker

    % Draw line from origin to the last valid point
    if ~isempty(validIndices)
        lastIdx = validIndices(end);
        plot3([0, pos(lastIdx, 1)], ...
              [0, pos(lastIdx, 2)], ...
              [0, pos(lastIdx, 3)], 'k-o');
    end

    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    hold off;
end