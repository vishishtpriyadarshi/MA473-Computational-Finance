function [U, S, t] = CrankNicolson(N, M, x, t, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, option)
    % This function solves the given parabolic PDE using CrankNicolson scheme.
    % INPUT: N, M, x, t, lambda, K, T, sigma, q, q_delta, g1, g2, phi
    % where, u(x, 0) = phi(x), u(a, t) = g1(t), u(b, t) = g2(t)
    % OUTPUT: U (representing the solution matrix)
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = phi(x, q_delta);
    U(:, 1) = g1(x(1), t, q_delta);
    U(:, end) = g2(x(end), t, q_delta);
    
    % Prepare tri-diagonal matrix
    A = tridiagonal_matrix(N - 1, 1 + 2*lambda, -lambda, -lambda);    
    L = tridiagonal_matrix(N - 1, 1 - 2*lambda, lambda, lambda);
    
    % Compute the values of U(x, t) at grid points
    for j = 2 : M + 1
        b = U(j-1, 2:end-1) + U(j-1, 1:end-2)*lambda/2 + U(j-1,2:end-1)*(1-lambda) + U(j-1,3:end)*lambda/2;         
        g_val = g(x(2 : end - 1), t(j), q, q_delta);
        U(j, 2 : end - 1) = psor_method(A, b - A * g_val', 1000, 1e-5)' + g_val;
    end
    
    % Get solution matrix for original eqn
    [U, S, t] = transform_pde(U, x, t, K, T, sigma, q, q_delta, option);
end