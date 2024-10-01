function printTraj(pos,posDot,N_start,N_end,color)
    figure;
    grid on, hold on
    if nargin > 4
        colormap(color);
    else
        colormap(winter);
    end
    speed = vecnorm(posDot')';
    patch([pos(N_start:N_end,1);nan],[pos(N_start:N_end,2);nan],[pos(N_start:N_end,3); nan],[speed(N_start:N_end); nan],'FaceColor','none','EdgeColor','interp','LineWidth', 1.5);
    cb = colorbar;
    ylabel(cb,'Ground Speed (m/s)','Rotation',270)
    view(3)
    plot3(0,0,0,'k*')
    plot3([0, pos(N_end,1)],...
         [0, pos(N_end,2)],...
         [0, pos(N_end,3)], 'k-o')
    xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)'), hold off
end