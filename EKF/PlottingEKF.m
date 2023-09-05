
Nmax = 1200;

%% Position
figure(1), title('Position')
subplot(3,1,1)
plot(Position.time(1:Nmax), xhatV(1:Nmax,1)), grid on, hold on
plot(Position.time(1:Nmax), Position.signals.values(1:Nmax,1))
legend('$r[1]$', '$\hat{r}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(Position.time(1:Nmax), xhatV(1:Nmax,2)), grid on, hold on
plot(Position.time(1:Nmax), Position.signals.values(1:Nmax,2))
legend('$r[2]$', '$\hat{r}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(Position.time(1:Nmax), xhatV(1:Nmax,3)), grid on, hold on
plot(Position.time(1:Nmax), Position.signals.values(1:Nmax,3))
legend('$r[3]$', '$\hat{r}[3]$', 'Interpreter', 'latex')


%% Speed
figure(2), title('Speed')
subplot(3,1,1)
plot(Position.time(1:Nmax), xhatV(1:Nmax,4)), grid on, hold on
plot(Position.time(1:Nmax), PositionDot.signals.values(1:Nmax,1))
legend('$\dot{r}[1]$', '$\dot{\hat{r}}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(Position.time(1:Nmax), xhatV(1:Nmax,5)), grid on, hold on
plot(Position.time(1:Nmax), PositionDot.signals.values(1:Nmax,2))
legend('$\dot{r}[2]$', '$\dot{\hat{r}}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(Position.time(1:Nmax), xhatV(1:Nmax,6)), grid on, hold on
plot(Position.time(1:Nmax), PositionDot.signals.values(1:Nmax,3))
legend('$\dot{r}[3]$', '$\dot{\hat{r}}[3]$', 'Interpreter', 'latex')

%% Acceleration
figure(3), title('Acceleration')
subplot(3,1,1)
plot(Position.time(1:Nmax), xhatV(1:Nmax,7)), grid on
legend('$\ddot{r}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(Position.time(1:Nmax), xhatV(1:Nmax,8)), grid on
legend('$\ddot{r}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(Position.time(1:Nmax), xhatV(1:Nmax,9)), grid on
legend('$\ddot{r}[3]$', 'Interpreter', 'latex')

%% Lagrange Multiplier
figure(4), title('Lagrange Multiplier')
plot(Position.time(1:Nmax), xhatV(1:Nmax,9)), grid on
legend('$\nu$', 'Interpreter', 'latex')

%% Nominal Wind


%% Total Force







