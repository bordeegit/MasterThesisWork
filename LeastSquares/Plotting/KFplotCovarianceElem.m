function KFplotCovarianceElem(P)

N = length(P);

figure;

% Extract matrix elements over time
p11 = squeeze(P(1,1,:));
p12 = squeeze(P(1,2,:));
p21 = squeeze(P(2,1,:));
p22 = squeeze(P(2,2,:));

% Time vector
t = 1:N;

% Plot all elements
subplot(2,2,1);
plot(t, p11, 'LineWidth', 1.5);
title('$P_{11}$ over time');
grid on;

subplot(2,2,2);
plot(t, p12, 'LineWidth', 1.5);
title('$P_{12}$ over time');
grid on;

subplot(2,2,3);
plot(t, p21, 'LineWidth', 1.5);
title('$P_{21}$ over time');
grid on;

subplot(2,2,4);
plot(t, p22, 'LineWidth', 1.5);
title('$P_{22}$ over time');
grid on;