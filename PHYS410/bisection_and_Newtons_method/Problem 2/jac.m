function A = jac(x0)

    x = sym('x', [1 length(x0)]);
    func = @(x) f(x);
    J = jacobian(func(x), x);
    
    A = eval(subs(J, x, x0));
   
end