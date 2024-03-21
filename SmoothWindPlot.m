
Nx = 10; % Number of points in X direction
Ny = 10; % Number of points in Y direction
Nz = 10; % Number of points in Z direction

[x_points, y_points, z_points] = meshgrid(linspace(0,100,Nx), linspace(-50,50,Ny), linspace(10,40,Nz));

% Creating wind vectors and populating (careful of limits)
WindVectors = zeros(Nx, Ny, Nz, 3);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            vx = WxFunctionZ(z_points(i,j,k));
            vy = WyFunctionXY(x_points(i,j,k), y_points(i,j,k));
            WindVectors(i,j,k,1) = vx; % x-component
            WindVectors(i,j,k,2) = vy; % y-component
        end
    end
end

% 3d Plot
quiver3(x_points, y_points, z_points, WindVectors(:,:,:,1), WindVectors(:,:,:,2), WindVectors(:,:,:,3), 'b');

% Wx

figure
plot(w_x_height_Data,height_Data,'o','LineWidth',1), hold on;
plot(WxFunctionZ(1:50),1:50, 'LineWidth', 1), hold off
legend('Data Points', 'Interpolated Curve', 'Location','southeast', 'Interpreter', 'latex');
ylabel('Height (m)', 'Interpreter', 'latex');
xlabel('Speed (m/s)', 'Interpreter', 'latex');
title('Wind in X direction', 'Interpreter', 'latex');

%ax = gca; 
%exportgraphics(ax, 'wxgen.pdf', 'ContentType','vector');

% Wy

x_pos_Data = [0,50,30,60,60,60,25,0,0,100,100,100,80];
y_pos_Data = [0,10,-10,0,-50,50,0,-50,50,0,-50,50,-20];
w_y_pos_Data = [0,10,5,0,0,0,0,0,0,0,0,0,-10];

[Xplot,Yplot] = meshgrid(0:2:120,-60:2:60);

figure
plot3(x_pos_Data,y_pos_Data,w_y_pos_Data, '.', 'MarkerSize', 30, 'Color', 'black'), hold on
surf(Xplot,Yplot, WyFunctionXY(Xplot, Yplot)), hold off
xlabel('X (m)', 'Interpreter', 'latex'), ylabel('Y (m)', 'Interpreter', 'latex');
zlabel('Speed (m/s)', 'Interpreter', 'latex');
title('Wind in Y direction' , 'Interpreter', 'latex'), grid on

%ax = gca; 
%exportgraphics(ax, 'wygen.pdf', 'ContentType','vector');


%% Contour 
for height = [10,15,20,25]
    % Specify the desired height
    desired_height = height; % in meters
    
    % Find the index of the closest point in the z direction to the desired height
    [~, z_index] = min(abs(z_points(1, 1, :) - desired_height));
    
    % Extract the wind vectors at the desired height
    vx_at_height = squeeze(WindVectors(:, :, z_index, 1));
    vy_at_height = squeeze(WindVectors(:, :, z_index, 2));
    
    % Calculate the magnitude of the wind at the desired height
    wind_magnitude_at_height = sqrt(vx_at_height.^2 + vy_at_height.^2);
    
    % Create a contour plot
    figure;
    contourf(x_points(:,:,1), y_points(:,:,1), wind_magnitude_at_height,10);
    clim([5 15]);
    colorbar;
    title(['Wind Magnitude at Height ', num2str(desired_height), 'm'], Interpreter="latex");
    xlabel('X (m)', Interpreter='latex');
    ylabel('Y (m)', Interpreter='latex');
    
    ax = gca; 
    exportgraphics(ax, append("cotour",num2str(desired_height),".png"));

end



