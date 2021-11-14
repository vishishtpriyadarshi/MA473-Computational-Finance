function x = gauss_seidel_method(A, b, iterations_count, epsilon)
    % This function solves the system of equation Ax = b using Gauss-Seidel Iteration method.
    
    n = length(A);
    D = diag(A);
    L = tril(A)- eye(n) .* D;
    U = triu(A)- eye(n) .* D;
    
    x = zeros(n, 1);
    x_new = x;
    for i = 1 : iterations_count
        for j = 1 : n
			x_new(j) = (b(j) - U(j, :) * x - L(j, :) * x_new) ./ D(j);
        end
        
        if norm(x_new - x) < epsilon
            break;
        end
        x = x_new;
    end
end