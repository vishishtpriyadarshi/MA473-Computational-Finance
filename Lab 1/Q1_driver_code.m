schemes = { @(N, M, x, t, r, dx, dt, f, g1, g2, phi) FTCS(N, M, x, t, r, dx, dt, f, g1, g2, phi); 
           @(N, M, x, t, r, dx, dt, f, g1, g2, phi) BTCS(N, M, x, t, r, dx, dt, f, g1, g2, phi);
           @(N, M, x, t, r, dx, dt, f, g1, g2, phi) CrankNicolson(N, M, x, t, r, dx, dt, f, g1, g2, phi) };

for i = 1 : length(schemes)
    % Take input variables for the PDE
    [N, M, a, b, T, f, g1, g2, phi] = input_pde();
    
    % Create grid and get mesh parameters
    [x, t, dx, dt, r] = create_grid(N, M, a, b, T);
    
    % Apply suitable scheme to get the solution matrix
    U = schemes{i}(N, M, x, t, r, dx, dt, f, g1, g2, phi);
    
    % Error computation for the numerical solution
    [N_list, error, log_error] = calculate_error(a, b, T, f, g1, g2, phi, schemes{i}, 1);
    
    % Tabulate and display the required plots
    tabulate(N_list, error, log_error, i);
    surface_plot(x, t, U, i);
end

function [N, M, a, b, T, f, g1, g2, phi] = input_pde()
    % This is a helper function to take the input values.
    
    N = 10;
    M = 2 * (N^2);
    a = 0;
    b = 1;
    T = 1;
    
    g1 = @(t) t * 0;
    g2 = @(t) t * 0;
    phi = @(x) x * 0;
    f = @(x, t) sin(2*pi*x) * sin(4*pi*t);
end