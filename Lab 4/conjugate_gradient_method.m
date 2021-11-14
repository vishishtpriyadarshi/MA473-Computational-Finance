function x = conjugate_gradient_method(A, b, iterations_count, epsilon)
    % This function solves the system of equation Ax = b using Conjugate Gradient method.
    
    n = length(A);
    [p_curr, r_curr, r_old] = deal(0, b, b);
    
    x_curr = zeros(n, 1);
    for i = 1 : iterations_count
        if i == 1
            p_new = r_curr;
        else
            p_new = r_curr + (norm(r_curr) ^ 2 / norm(r_old) ^ 2) * p_curr;
        end
        
        alpha = norm(r_curr) ^ 2 / (p_new' * A * p_new);
        x_new = x_curr + alpha * p_new;
        r_new = r_curr - alpha * A * p_new;
        
        [r_old, r_curr, x_curr] = deal(r_curr, r_new, x_new);
        
        if norm(r_new) < epsilon
            break;
        end
    end
    
    x = x_curr;
end

