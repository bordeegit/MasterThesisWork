close all
Nmax = 5000;

timePlot = Position.time(1:Nmax);

%% Position
figure(1), sgtitle('Position')
subplot(3,1,1)
plot(timePlot, Position.signals.values(1:Nmax,1)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,1))
legend('$r[1]$', '$\hat{r}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(timePlot, Position.signals.values(1:Nmax,2)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,2))
legend('$r[2]$', '$\hat{r}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(timePlot, Position.signals.values(1:Nmax,3)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,3))
legend('$r[3]$', '$\hat{r}[3]$', 'Interpreter', 'latex')



%% Speed
figure(2), sgtitle('Speed')
subplot(3,1,1)
plot(timePlot, PositionDot.signals.values(1:Nmax,1)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,4))
legend('$\dot{r}[1]$', '$\dot{\hat{r}}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(timePlot, PositionDot.signals.values(1:Nmax,2)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,5))
legend('$\dot{r}[2]$', '$\dot{\hat{r}}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(timePlot, PositionDot.signals.values(1:Nmax,3)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,6))
legend('$\dot{r}[3]$', '$\dot{\hat{r}}[3]$', 'Interpreter', 'latex')

%% Acceleration
figure(3), sgtitle('Acceleration')
subplot(3,1,1)
plot(timePlot, PositionDotDot.signals.values(1:Nmax,1)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,7))
legend('$\ddot{r}[1]$', '$\ddot{\hat{r}}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(timePlot, PositionDotDot.signals.values(1:Nmax,2)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,8))
legend('$\ddot{r}[2]$', '$\ddot{\hat{r}}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(timePlot, PositionDotDot.signals.values(1:Nmax,3)), grid on, hold on
plot(timePlot, xhatV(1:Nmax,9))
legend('$\ddot{r}[3]$', '$\ddot{\hat{r}}[3]$', 'Interpreter', 'latex')

%% Nominal Wind
figure(4), sgtitle('Nominal Wind')
subplot(2,1,1)
plot(W_log.time(1:Nmax), W_log.signals.values(1:Nmax,1)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,10))
legend('$w_{nx}$', '$\hat{w_{nx}}$', 'Interpreter', 'latex')

subplot(2,1,2)
plot(W_log.time(1:Nmax), W_log.signals.values(1:Nmax,2)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,11))
legend('$w_{ny}$', '$\hat{w_{ny}}$', 'Interpreter', 'latex')

%% Lift Force

figure(5), sgtitle('Lift Force')
subplot(3,1,1)
plot(Forces.time(1:Nmax), Forces.signals.values(1:Nmax,1)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,12))
legend('$F_{totx}$', '$\hat{F}_{totx}$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(Forces.time(1:Nmax), Forces.signals.values(1:Nmax,2)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,13))
legend('$F_{toty}$', '$\hat{F}_{toty}$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(Forces.time(1:Nmax), Forces.signals.values(1:Nmax,3)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,14))
legend('$F_{totz}$', '$\hat{F}_{totz}$', 'Interpreter', 'latex')

%% Drag Force Magnitude

figure(6), sgtitle('Drag Force Magnitude')
plot(W_log.time(1:Nmax), xhatV(1:Nmax,15)), grid on
legend('$|F_{D}|$', 'Interpreter', 'latex')

%% Linear relationship between the steering input and the angular rate

figure(7), sgtitle('Steering Coefficient')
plot(W_log.time(1:Nmax), xhatV(1:Nmax,16)), grid on
legend('$c_u$', 'Interpreter', 'latex')



