close all

% Line force
figure(1),plot(Forces.time,Forces.signals.values(:,end)/1000,'k'),grid on,hold on
xlabel('time (s)'), ylabel('Pulling force (kN)')

% Path in spherical coord. force
figure(2),plot(phi.signals.values,theta.signals.values,'k'),grid on,hold on
xlabel('\phi (rad)'), ylabel('\theta (rad)')

% Wing speed
figure(3),plot(Vkite.time,Vkite.signals.values,'k'),grid on,hold on
xlabel('time (s)'), ylabel('Kite speed (m/s)')

% Angle of attack
figure(4),plot(alpha.time,alpha.signals.values*180/pi,'k'),grid on,hold on
xlabel('time (s)'), ylabel('Angle of attack (deg)')

% 3D path
figure(5),plot3(Position.signals.values(:,1),Position.signals.values(:,2),Position.signals.values(:,3),'k'),grid on,hold on
plot3(0,0,0,'k*')
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
WindPlot

% Theta, phi
figure(6),plot(theta.time,theta.signals.values(:,end),'k',phi.time,phi.signals.values(:,end),'k-.'),grid on,hold on
xlabel('time (s)'), ylabel('\theta,\phi (rad)'), legend('\theta','\phi')

% r
figure(7),plot(r.time,r.signals.values(:,end),'k'),grid on,hold on
xlabel('time (s)'), ylabel('r (m)')

% Wind
figure(8)
subplot(3,1,1)
plot(W_log.time, W_log.signals.values(:,1),'k'), grid on, hold on
xlabel('time (s)'), ylabel('Wind X (m/s)')
subplot(3,1,2)
plot(W_log.time, W_log.signals.values(:,2),'k'), grid on, hold on
xlabel('time (s)'), ylabel('Wind Y (m/s)')
subplot(3,1,3)
plot(W_log.time, W_log.signals.values(:,3),'k'), grid on, hold on
xlabel('time (s)'), ylabel('Wind Z (m/s)')


