y0 = [1; 1];

dt = 0.2;
y_exact = [exp(0.5*dt); exp(0.25*dt)];
yout = rk4step(@fcn, 0, dt, y0);
err1 = abs(y_exact - yout);

dt = 0.1;
y_exact = [exp(0.5*dt); exp(0.25*dt)];
yout = rk4step(@fcn, 0, dt, y0);
err2 = abs(y_exact - yout);

R = err1./err2;