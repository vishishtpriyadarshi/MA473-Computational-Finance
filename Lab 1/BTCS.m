function U = BTCS(N, M, x, t, r, dx, dt, f, g1, g2, phi)
    % This function solves the given parabolic PDE using BTCS scheme.
    % INPUT: N, M, x, t, r, dx, dt, f, g1, g2, phi
    % where, u(x, 0) = phi(x), u(0, t) = g1(t), u(x, t) = g2(t)
    % and the non-homogeneous function is f.
    % OUTPUT: U (representing the solution matrix)
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = phi(x);
    U(:, 1) = g1(t);
    U(:, end) = g2(t);
    
    % Prepare tri-diagonal matrix
    L = tridiagonal_matrix(N - 1, 1 + 2*r, -r, -r);
    
    % Compute the values of U(x, t) at grid points
    for j = 2 : M
        u = L \ (U(j - 1, 2 : end - 1)' + dt * f(x(2: end - 1), t(j))');
        U(j, 2 : end - 1) = u;
    end
end