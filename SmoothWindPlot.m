
Nx = 10; % Number of points in X direction
Ny = 10; % Number of points in Y direction
Nz = 10; % Number of points in Z direction

[x_points, y_points, z_points] = meshgrid(linspace(40,100,Nx), linspace(-50,50,Ny), linspace(10,40,Nz));

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

% Plot
quiver3(x_points, y_points, z_points, WindVectors(:,:,:,1), WindVectors(:,:,:,2), WindVectors(:,:,:,3), 'b');