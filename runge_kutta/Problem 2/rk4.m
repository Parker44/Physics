function [tout yout] = rk4(fcn, tspan, y0)
% Inputs
% fcn: Function handle for right hand sides of ODEs (returns
% length-n column vector)
% tspan: Vector of output times (length nout).
% y0: Initial values (length-n column vector).
%
% Outputs
% tout: Vector of output times.
% yout: Output values (nout x n array).
    
    tout = tspan;
    dt = tout(2) - tout(1);
    
    n = length(y0);
    nout = length(tspan);
    yout = zeros(nout, n);
    
    yout(1,:) = y0;
    
    for i=2:nout
        f0 = fcn(yout(i-1,:));
        f1 = fcn(yout(i-1,:) + (dt/2)*f0');
        f2 = fcn(yout(i-1,:) + (dt/2)*f1');
        f3 = fcn(yout(i-1,:) + dt*f2');

        yout(i,:) = yout(i-1,:) + (dt/6)*(f0' + 2*f1' + 2*f2' + f3');
    end
end

