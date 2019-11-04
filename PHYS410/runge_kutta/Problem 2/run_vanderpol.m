y0 = [1 -6];
tmax = 100;

clf;
hold on

level = 12;
nt = 2^level + 1;
tspan = linspace(0.0, tmax, nt);
[tout, yout] = rk4(@vanderpol, tspan, y0);

%plot(tout, yout(:, 1), 'r-');

plot(yout(:,1), yout(:, 2), 'b-');
