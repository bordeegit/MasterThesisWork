f = figure;
grid on, hold on
colormap(winter);
speed = vecnorm(PositionDot.signals.values')';
patch([Position.signals.values(:,1);nan],[Position.signals.values(:,2);nan],[Position.signals.values(:,3); nan],[speed; nan],'FaceColor','none','EdgeColor','interp','LineWidth', 1.5);
cb = colorbar;
ylabel(cb,'Ground Speed (m/s)','Rotation',270)
view(3)
plot3(0,0,0,'k*')
plot3([0, Position.signals.values(end,1)],...
     [0, Position.signals.values(end,2)],...
     [0, Position.signals.values(end,3)], 'k-o')
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)'), hold off
