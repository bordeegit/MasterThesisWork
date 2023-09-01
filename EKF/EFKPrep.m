%% PositionDot computation
% Creation of the (paper) rdot vector, that in this case is called 
% PositionDot since 'rdot' is the name of Ldot in the paper.
% Equations are derived in (16) of the paper.
% PositionDot is created based on Position, so that they will share the
% same structure, as if it was logged during the simulation.
Plotting = false;
PositionDot = Position;
PositionDot.signals.label = '<PositionDot>';

L = r.signals.values;
Ldot = rdot.signals.values;
th = theta.signals.values;
thdot = thetadot.signals.values;
ph = phi.signals.values;
phdot = phi.signals.values;

PositionDot.signals.values = [Ldot.*sin(th).*cos(ph)+L.*(thdot.*cos(th).*cos(ph)-phdot.*sin(th).*sin(ph)) ... 
               Ldot.*sin(th).*sin(ph)+L.*(thdot.*cos(th).*sin(ph)+phdot.*sin(th).*cos(ph)) ...
               Ldot.*cos(th)-L.*thdot.*sin(th)];

% Plotting
if(Plotting)
figure
subplot(3,1,1)
plot(PositionDot.time, PositionDot.signals.values(:,1),'k'), grid on
xlabel('time (s)'), ylabel('Kite Vel X (m/s)')
subplot(3,1,2)
plot(PositionDot.time, PositionDot.signals.values(:,2),'k'), grid on
xlabel('time (s)'), ylabel('Kite Vel Y (m/s)')
subplot(3,1,3)
plot(PositionDot.time, PositionDot.signals.values(:,3),'k'), grid on
xlabel('time (s)'), ylabel('Kite Vel Z (m/s)')
end
