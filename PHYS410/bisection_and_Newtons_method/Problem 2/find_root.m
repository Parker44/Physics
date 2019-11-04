x = [3 -2 -1];
tol = 10^(-16);

x = newtond(@f, @jac, x, tol);