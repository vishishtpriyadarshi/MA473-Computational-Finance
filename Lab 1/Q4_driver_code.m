schemes = { @(N, M, x, t, r, dx, dt, f, g1, g2, phi) FTCS_modified(N, M, x, t, r, dx, dt, f, g1, g2, phi); 
           @(N, M, x, t, r, dx, dt, f, g1, g2, phi) BTCS_modified(N, M, x, t, r, dx, dt, f, g1, g2, phi);
           @(N, M, x, t, r, dx, dt, f, g1, g2, phi) CrankNicolson_modified(N, M, x, t, r, dx, dt, f, g1, g2, phi) };

for i = 1 : length(schemes)
    % Take input variables for the PDE
    [N, M, a, b, T, f, g1, g2, phi] = input_pde();
    
    % Create grid and get mesh parameters
    [x, t, dx, dt, r] = create_grid(N, M, a, b, T);
    
    % Apply suitable scheme to get the solution matrix
    U = schemes{i}(N, M, x, t, r, dx, dt, f, g1, g2, phi);
    
    % Error computation for the numerical solution
    [N_list, error, log_error] = calculate_error(a, b, T, f, g1, g2, phi, schemes{i}, 2);
    
    % Tabulate and display the required plots
    tabulate(N_list, error, log_error, i);
    surface_plot(x, t, U, i);
end

function [N, M, a, b, T, f, g1, g2, phi] = input_pde()
    % This is a helper function to take the input values.
    
    N = 10;
    M = 2 * (N^2);
    a = 0;
    b = 1/2;
    T = 1;
    
    g1 = @(t) t * 0;
    g2 = @(t) t .* t;
    phi = @(x) x .* (1 - x);
    f = @(x, t) x * 0;
end

function U = FTCS_modified(N, M, x, t, r, dx, dt, f, g1, g2, phi)
    % Modified FTCS scheme covering Neumann conditions
    
    if r > 0.5
        error("Scheme is unstable - r > 0.5");
    end
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = phi(x);
    U(:, 1) = g1(t);
    
    % Prepare tri-diagonal matrix
    L = tridiagonal_matrix(N - 1, 1 - 2*r, r, r);
    L(end, end) = 1 - r;
    b = zeros(N - 1, 1);
    
    % Compute the values of U(x, t) at grid points
    for j = 2 : M + 1
        b(end) = r * dx * g2(t(j - 1));
        u = L * U(j - 1, 2 : end - 1)' + b + dt * f(x(2: end - 1), t(j - 1))';
        U(j, 2 : end - 1) = u;
    end
    U(:, end) = U(:, end - 1) + (dx * g2(t))';
end

function U = BTCS_modified(N, M, x, t, r, dx, dt, f, g1, g2, phi)
    % Modified BTCS scheme covering Neumann conditions
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = phi(x);
    U(:, 1) = g1(t);
    
    % Prepare tri-diagonal matrix
    L = tridiagonal_matrix(N - 1, 1 + 2*r, -r, -r);
    L(end, end) = 1 + r;
    b = zeros(N - 1, 1);
    
    % Compute the values of U(x, t) at grid points
    for j = 2 : M + 1
        b(end) = -r * dx * g2(t(j - 1));
        u = L \ (U(j - 1, 2 : end - 1)' - b + dt * f(x(2: end - 1), t(j))');
        U(j, 2 : end - 1) = u;
    end
    U(:, end) = U(:, end - 1) + (dx * g2(t))';
end

function U = CrankNicolson_modified(N, M, x, t, r, dx, dt, f, g1, g2, phi)
    % Modified Crank-Nicolson scheme covering Neumann conditions
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    U(1, :) = phi(x);
    U(:, 1) = g1(t);
    
    % Prepare tri-diagonal matrix
    A = tridiagonal_matrix(N - 1, 1 + 2*r, -r, -r);  
    A(end, end) = 1 + r;    
    L = tridiagonal_matrix(N - 1, 1 - 2*r, r, r);
    L(end, end) = 1 - r;
    
    b1 = zeros(N - 1, 1);
    b2 = zeros(N - 1, 1);
    
    % Compute the values of U(x, t) at grid points
    for j = 2 : M + 1
        b1(end) = r * dx * g2(t(j - 1));
        b2(end) = -r * dx * g2(t(j));
        u = A \ (L*(U(j - 1, 2 : end - 1)' + b1 - b2 + dt * f(x(2: end - 1), t(j - 1))'));
        U(j, 2 : end - 1) = u;
    end
    U(:, end) = U(:, end - 1) + (dx * g2(t))';
end