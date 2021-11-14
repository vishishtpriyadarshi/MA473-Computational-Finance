function U = FTCS(N, M, x, t, lambda, dx, dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m)
    % This function solves the given general parabolic PDE using FTCS scheme.
    % OUTPUT: U (representing the solution matrix)
 
    if lambda > 0.5
        error("Scheme is unstable - lambda > 0.5");
    end
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = g(x, K);
    U(:, 1) = f1(x(1), t, K, r, T);
    U(:, end) = f2(x(end), t, K, r, T);
    
    B = zeros(N - 1, 1);
    
    % Prepare tri-diagonal matrix
    L = tridiagonal_matrix(N - 1, x, t, r, sigma, delta, b, c, d, dx, dt, 1);
        
    % Compute the values of U(x, t) at grid points
    for j = 2 : M + 1
        alpha = -b(x(1), t, sigma) * dt / dx ^ 2 + 0.5 * c(x(1), t, r, delta) * dt / dx;
        gamma = -b(x(end - 1), t, sigma) * dt / dx ^ 2 - 0.5 * c(x(end - 1), t, r, delta) * dt / dx;
        B(1) = alpha * U(j - 1, 1); 
        B(end) = gamma * U(j - 1, end);
        
        u = L * (U(j - 1, 2 : end - 1)') + B;
        U(j, 2 : end - 1) = u;
    end

    U = flipud(U);
end