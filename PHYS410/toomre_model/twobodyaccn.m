function [a] = twobodyaccn(m1, m2, r)
    
    a = zeros(2, 3);
    
    deltar = r(2, :) - r(1, :);
    a(1, :) = (m2/norm(deltar)^3)*deltar;
    
    deltar = r(1, :) - r(2, :);
    a(2, :) = (m1/norm(deltar)^3)*deltar;
    
end

