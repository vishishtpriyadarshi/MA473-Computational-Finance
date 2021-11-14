function x = psor_method(A, b, iterations_count, epsilon)
    % This function solves the system of equation Ax = b using PSOR Iteration method.
    
    n = length(A);
    D = diag(A);
    L = tril(A)- eye(n) .* D;
    U = triu(A)- eye(n) .* D;
    
    Bj = -1 * (eye(n) ./ D) * (L + U);
    w = 2 / (1 + sqrt(1 - norm(Bj) ^ 2));
    
    x = zeros(n, 1);
    x_new = x;
    for i = 1 : iterations_count
        for j = 1 : n
			x_new(j) = max(0, (1 - w) * x(j) + w * (b(j) - U(j, :) * x - L(j, :) * x_new) ./ D(j));
        end
        
        if norm(x_new - x) < epsilon
            break;
        end
        x = x_new;
    end
end

