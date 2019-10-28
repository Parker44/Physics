function m = dfdx(x0)
    
    syms x
    func = @(x) f(x);
    m = eval(subs(diff(func, x), x, x0));
    
end