y0 = [0 1];
tmax = 3*pi;

clf;
hold on

level = 6;
nt = 2^level + 1;
tspan = linspace(0.0, tmax, nt);
[tout6, yout6] = rk4(@shm, tspan, y0);


level = 7;
nt = 2^level + 1;
tspan = linspace(0.0, tmax, nt);
[tout7, yout7] = rk4(@shm, tspan, y0);


level = 8;
nt = 2^level + 1;
tspan = linspace(0.0, tmax, nt);
[tout8, yout8] = rk4(@shm, tspan, y0);

level = 9;
nt = 2^level + 1;
tspan = linspace(0.0, tmax, nt);
[tout9, yout9] = rk4(@shm, tspan, y0);

yout7 = yout7(1:2:end, :);
err67 = yout6-yout7;

yout8 = yout8(1:4:end, :);
err78 = yout7-yout8;

yout9 = yout9(1:8:end, :);
err89 = yout8-yout9;

plot(tout6, err67(:,1), 'r-');
plot(tout6, (2^4)*err78(:,1), 'g-');
plot(tout6, (4^4)*err89(:,1), 'b-');