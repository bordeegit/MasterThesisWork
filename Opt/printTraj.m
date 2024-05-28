figure;
grid on, hold on
colormap(winter);
speed = vecnorm(posDot')';
patch([pos(:,1);nan],[pos(:,2);nan],[pos(:,3); nan],[speed; nan],'FaceColor','none','EdgeColor','interp','LineWidth', 1.5);
cb = colorbar;
ylabel(cb,'Ground Speed (m/s)','Rotation',270)
view(3)
plot3(0,0,0,'k*')
plot3([0, pos(end,1)],...
     [0, pos(end,2)],...
     [0, pos(end,3)], 'k-o')
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)'), hold off
