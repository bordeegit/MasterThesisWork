%% PositionDot computation
% Creation of the (paper) rdot vector, that in this case is called 
% PositionDot since 'rdot' is the name of Ldot in the paper.
% Equations are derived in (16) of the paper.
% PositionDot is created based on Position, so that they will share the
% same structure, as if it was logged during the simulation.

% Note that theta and thetadot are converted back to the convention used in
% the model (they where modified because of the controller)

Plotting = false;
PositionDot = Position;
PositionDot.signals.label = '<PositionDot>';

L = r.signals.values;
Ldot = rdot.signals.values;
th = pi/2 - theta.signals.values;  % Reverted back
thdot = -thetadot.signals.values;  % Reverted back
ph = phi.signals.values;
phdot = phi.signals.values;

PositionDot.signals.values = [Ldot.*sin(th).*cos(ph)+L.*(thdot.*cos(th).*cos(ph)-phdot.*sin(th).*sin(ph)) ... 
               Ldot.*sin(th).*sin(ph)+L.*(thdot.*cos(th).*sin(ph)+phdot.*sin(th).*cos(ph)) ...
               Ldot.*cos(th)-L.*thdot.*sin(th)];

%% Computation of drag magnitude 

sz = size(Forces.time,1);

drag_norm = zeros(sz,1);

for i = 1:sz
    drag_norm(i) = norm(Forces.signals.values(i,4:6));
end


%% Plotting
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

