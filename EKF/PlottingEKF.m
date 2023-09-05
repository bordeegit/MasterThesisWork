
Nmax = 1200;

%% Position
figure(1), title('Position')
subplot(3,1,1)
plot(Position.time(1:Nmax), Position.signals.values(1:Nmax,1)), grid on, hold on
plot(Position.time(1:Nmax), xhatV(1:Nmax,1))
legend('$r[1]$', '$\hat{r}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(Position.time(1:Nmax), Position.signals.values(1:Nmax,2)), grid on, hold on
plot(Position.time(1:Nmax), xhatV(1:Nmax,2))
legend('$r[2]$', '$\hat{r}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(Position.time(1:Nmax), Position.signals.values(1:Nmax,3)), grid on, hold on
plot(Position.time(1:Nmax), xhatV(1:Nmax,3))
legend('$r[3]$', '$\hat{r}[3]$', 'Interpreter', 'latex')


%% Speed
figure(2), title('Speed')
subplot(3,1,1)
plot(Position.time(1:Nmax), PositionDot.signals.values(1:Nmax,1)), grid on, hold on
plot(Position.time(1:Nmax), xhatV(1:Nmax,4))
legend('$\dot{r}[1]$', '$\dot{\hat{r}}[1]$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(Position.time(1:Nmax), PositionDot.signals.values(1:Nmax,2)), grid on, hold on
plot(Position.time(1:Nmax), xhatV(1:Nmax,5))
legend('$\dot{r}[2]$', '$\dot{\hat{r}}[2]$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(Position.time(1:Nmax), PositionDot.signals.values(1:Nmax,3)), grid on, hold on
plot(Position.time(1:Nmax), xhatV(1:Nmax,6))
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
plot(Position.time(1:Nmax), xhatV(1:Nmax,10)), grid on
legend('$\nu$', 'Interpreter', 'latex')

%% Nominal Wind
figure(5), title('Nominal Wind')
subplot(2,1,1)
plot(W_log.time(1:Nmax), W_log.signals.values(1:Nmax,1) - WindX), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,11))
legend('$w_{nx}$', '$\hat{w_{nx}}$', 'Interpreter', 'latex')

subplot(2,1,2)
plot(W_log.time(1:Nmax), W_log.signals.values(1:Nmax,2)- WindY), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,12))
legend('$w_{ny}$', '$\hat{w_{ny}}$', 'Interpreter', 'latex')

%% Total Force
%F_LDg = Forces.signals.values(1:Nmax,1) + Forces.signals.values(1:Nmax,1);

figure(6), title('Total Force')
subplot(3,1,1)
plot(Forces.time(1:Nmax), Forces.signals.values(1:Nmax,1)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,13))
legend('$w_{nx}$', '$\hat{w}_{nx}$', 'Interpreter', 'latex')

subplot(3,1,2)
plot(Forces.time(1:Nmax), Forces.signals.values(1:Nmax,2)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,14))
legend('$w_{ny}$', '$\hat{w}_{ny}$', 'Interpreter', 'latex')

subplot(3,1,3)
plot(Forces.time(1:Nmax), Forces.signals.values(1:Nmax,3)), grid on, hold on
plot(W_log.time(1:Nmax), xhatV(1:Nmax,15))
legend('$w_{nz}$', '$\hat{w}_{nz}$', 'Interpreter', 'latex')






