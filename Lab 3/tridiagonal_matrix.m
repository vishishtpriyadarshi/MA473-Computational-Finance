function A = tridiagonal_matrix(N, x, t, r, sigma, delta, b, c, d, dx, dt, sgn)
    % This function creates a NxN tri-diagonal matrix.
    
    A = zeros(N);
    
    for i = 2 : N + 1
        if i > 2
            A(i - 1, i - 2) = sgn * -b(x(i), t, sigma) * dt / dx ^ 2 + sgn * 0.5 * c(x(i), t, r, delta) * dt / dx;
        end

        A(i - 1, i - 1) = 1 + sgn * 2 * b(x(i), t, sigma) * dt / dx ^ 2 - sgn * d(x(i), t, r) * dt;
        
        if i <= N
            A(i - 1, i) = sgn * -b(x(i), t, sigma) * dt / dx ^ 2 - sgn * 0.5 * c(x(i), t, r, delta) * dt / dx;
        end
    end
end

