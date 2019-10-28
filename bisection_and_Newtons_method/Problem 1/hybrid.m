function x = hybrid(f, dfdx, xmin, xmax, tol1, tol2)
    
    fmin = f(xmin);
    
    converged1 = 0;
    while ~converged1
        xmid = (xmin + xmax)/2;
        fmid = f(xmid);
        
        if fmid == 0
            break;
        elseif fmid*fmin < 0
            xmax = xmid;
        else
            xmin = xmid;
            fmin = fmid;
        end
        
        if (xmax - xmin)/abs(xmid) < tol1
            converged1 = 1;
        end
    end
    
    prev = xmid;

    converged2 = 0;
    while ~converged2
        dx = f(prev)/dfdx(prev);
        next = prev - dx;
        
        if dx/next < tol2
            converged2 = 1;
        else
            prev = next;
        end
    end
    
    x = next;
    
end