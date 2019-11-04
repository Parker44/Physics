function [tout, yout] = rk4ad(fcn, tspan, reltol, y0)
% Inputs
% fcn: Function handle for right hand sides of ODEs (returns
% length-n column vector)
% tspan: Vector of output times (length nout vector).
% reltol: Relative tolerance parameter.
% y0: Initial values (length-n column vector).
%
% Outputs
% tout: Output times (length-nout column vector, elements
% identical to tspan).
% yout: Output values (nout x n array).

    n = length(y0);
    tout = tspan;
    
    dtspan = tspan(2) - tspan(1);
    nout = length(tspan);
    yout = zeros(nout, n);
    
    yout(1,:) = y0;
    
    % For keeping track of current time during integration
    t = tspan(1);
    for i=2:nout
        N = 1;  % time step dividing factor to be incremented
        y_tmp = yout(i-1, :); % current y0 for rk4step
        while t < tspan(i)
            dt = tspan(i) - t;
            y_C = rk4step(fcn, 0, dt, y_tmp);
            
            y_F = rk4step(fcn, 0, dt/2, y_tmp);
            y_F = rk4step(fcn, 0, dt/2, y_F);
            
            e_C = (16/15)*abs(y_C(1) - y_F(1));
            
            while e_C > reltol
                N = N+1;
                dt = dtspan/N;
                if dt <= 1.0e-4
                    break;  % Reached minimum allowed time step
                end
                
                y_C = rk4step(fcn, 0, dt, y_tmp);
            
                y_F = rk4step(fcn, 0, dt/2, y_tmp);
                y_F = rk4step(fcn, 0, dt/2, y_F);

                e_C = (16/15)*abs(y_C(1) - y_F(1));
            end
            y_tmp = y_C;
            t = t + dt;
        end
        yout(i, :) = y_tmp;
    end

end