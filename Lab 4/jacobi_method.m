function x = jacobi_method(A, b, iterations_count, epsilon)
    % This function solves the system of equation Ax = b using Jacobi Iteration method.
    
    n = length(A);
    D = eye(n) .* diag(A);
    
    x = zeros(n, 1);
    for i = 1 : iterations_count
        x_new = (b - (A - D) * x) ./ diag(A);
        if norm(x_new - x) < epsilon
            break;
        end
        x = x_new;
    end
end

