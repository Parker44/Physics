format long

tol1 = 10^(-2);
tol2 = 10^(-16);

roots = zeros(1, 8);

xmin = -1;
dx = 0.1;

for i=1:length(roots)
    newroot = 0;
    while ~newroot
        xmax = xmin + dx;
        if f(xmin)*f(xmax) < 0
            roots(i) = hybrid(@f, @dfdx, xmin, xmax, tol1, tol2);
            newroot = 1;
        end
        xmin = xmax;
    end
end

