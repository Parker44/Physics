function [a] = nbodyaccn(cores_m, n_cores, cores_ns, r)
    a = zeros(size(r(:, :, :, 1)));
    for i=1:n_cores
        % Only the cores contribute to the gravitational field
        for j=1:cores_ns+1
            for k=1:n_cores
                if ~(i==k && j==1)
                    deltar = r(k,1,:) - r(i,j,:);
                    deltar = reshape(deltar, [1,3]);
                    a(i,j,:) = reshape(a(i,j,:), [1,3]) + (cores_m(k)/norm(deltar)^3)*deltar;
                end
            end
        end
    end
end