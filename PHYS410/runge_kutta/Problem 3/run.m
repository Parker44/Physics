% Perform rk4ad on shm function

clf;
hold on

tspan = linspace(0.0, 3.0 * pi, 65);
reltols = [1.0e-5 1.0e-7 1.0e-9 1.0e-11];
y0 = [0 1];

y_exact = sin(tspan');

[tout, yout] = rk4ad(@shm, tspan, reltols(4), y0);
err1 = yout(:, 1) - y_exact;

%plot(tout, err1, 'r-');

% Perform rk4ad on Van der Pol
tspan = linspace(0.0, 100, 4097);
reltol = 1.0e-10;

y0 = [1 -6];

[tout, yout] = rk4ad(@vanderpol, tspan, reltol, y0);

plot(tout, yout(:, 1), 'r-');
