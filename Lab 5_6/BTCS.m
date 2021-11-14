function [U, S, t] = BTCS(N, M, x, t, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, option)
    % This function solves the given parabolic PDE using BTCS scheme.
    % INPUT: N, M, x, t, lambda, K, T, sigma, q, q_delta, g1, g2, phi
    % where, u(x, 0) = phi(x), u(a, t) = g1(t), u(b, t) = g2(t)
    % OUTPUT: U (representing the solution matrix)
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = phi(x, q_delta);
    U(:, 1) = g1(x(1), t, q_delta);
    U(:, end) = g2(x(end), t, q_delta);
    b = zeros(N - 1, 1);
    
    % Prepare tri-diagonal matrix
    L = tridiagonal_matrix(N - 1, 1 + 2*lambda, -lambda, -lambda);
    
    % Compute the values of U(x, t) at grid points
    for j = 2 : M + 1
        if strcmp(option, 'call')
            b = U(j-1, 2 : end-1);
        end
        g_val = g(x(2 : end - 1), t(j), q, q_delta);
        U(j, 2 : end - 1) = psor_method(L, b - L * g_val', 1000, 1e-5)' + g_val;
    end
    
    % Get solution matrix for original eqn
    [U, S, t] = transform_pde(U, x, t, K, T, sigma, q, q_delta, option);
end