function yout = rk4step(fcn, ~, dt, y0)
% Inputs
% fcn: Function handle for right hand sides of ODEs (returns
% length-n column vector).
% t0: Initial value of independent variable.
% dt: Time step.
% y0: Initial values (length-n column vector).
%
% Output
% yout: Final values (length-n column vector).

    f0 = fcn(y0);
    f1 = fcn(y0 + (dt/2)*f0);
    f2 = fcn(y0 + (dt/2)*f1);
    f3 = fcn(y0 + dt*f2);

    yout = y0 + (dt/6)*(f0 + 2*f1 + 2*f2 + f3);
end

