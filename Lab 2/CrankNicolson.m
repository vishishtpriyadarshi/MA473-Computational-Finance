function [U, S, t] = CrankNicolson(N, M, x, t, lambda, K, T, sigma, q, q_delta, g1, g2, phi)
    % This function solves the given parabolic PDE using CrankNicolson scheme.
    % INPUT: N, M, x, t, lambda, K, T, sigma, q, q_delta, g1, g2, phi
    % where, u(x, 0) = phi(x), u(a, t) = g1(t), u(b, t) = g2(t)
    % OUTPUT: U (representing the solution matrix)
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = phi(x);
    U(:, 1) = g1(x(1), t);
    U(:, end) = g2(x(end), t);
    
    b = zeros(N - 1, 1);
    
    % Prepare tri-diagonal matrix
    A = tridiagonal_matrix(N - 1, 1 + 2*lambda, -lambda, -lambda);    
    L = tridiagonal_matrix(N - 1, 1 - 2*lambda, lambda, lambda);
    
    % Compute the values of U(x, t) at grid points
    for j = 2 : M + 1
        b(1) = lambda * (g1(x(1), t(j - 1)) + g1(x(1), t(j))); 
        b(end) = lambda * (g2(x(end), t(j - 1)) + g2(x(end), t(j)));
        
        u = A \ (L * (U(j - 1, 2 : end - 1)') + b);
        U(j, 2 : end - 1) = u;
    end
    
    % Get solution matrix for original eqn
    [U, S, t] = transform_pde(U, x, t, K, T, sigma, q, q_delta);
end


