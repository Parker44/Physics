function y = f(x)
    
    y = [x(1)^2 + x(2)^3 + x(3)^4 - 1;
         sin(x(1)*x(2)*x(3)) - x(1) - x(2) - x(3);
         x(1) - x(2)*x(3)];
    
end