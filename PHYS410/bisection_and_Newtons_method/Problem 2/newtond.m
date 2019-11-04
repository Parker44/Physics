function x = newtond(f, jac, x0, tol)

    x = x0;
    
    for a=1:100
       res = f(x);
       
       dx = (jac(x)\res)';
       x = x-dx;
       
       if norm(dx.^2)/sqrt(length(dx)) <= tol
           break;
       end
    end

end