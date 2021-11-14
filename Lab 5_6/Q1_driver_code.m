schemes = { @(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, option) BTCS(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, option);
            @(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, option) CrankNicolson(N, M, x, tau, lambda, K, T, sigma, q, q_delta,g,  g1, g2, phi, option)}; 

scheme_names = ["BTCS", "Crank-Nicolson"];

for i = 1 : length(schemes)
    % Take input variables for the PDE
    [h, k, a, b, K, T, sigma, q, q_delta, r, g, g1, g2, phi] = input_pde();
    
    % Create grid and get mesh parameters
    [x, tau, lambda, N, M] = create_grid(h*i, k*i, a, b, T * (sigma ^ 2) / 2);
    fprintf("\n\nExecuting %s scheme for Surface plot\n\n", scheme_names(i));
    
    % part(a): Apply suitable scheme to get the solution matrix
    [V, S, t] = schemes{1}(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, 'put');
    
    % part(b): Display the surface plot
    surface_plot(S, t, V, i, K, 'put');
    
    % part(c): Plot the error vs S plot
    [x, tau, lambda, N, M] = create_grid(h*i/2, k*i/2, a, b, T * (sigma ^ 2) / 2);
    [V2, S2, t2] = schemes{1}(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, 'put');
    V2 = V2(1:2:end, 1:2:end);
    
    figure();
    plot(S, abs(V(end, :) - V2(end, :)));
    xlabel("S");
    ylabel("Error");
    title("Error in V(S, t) using δx, δt and δx/2, δt/2: " + scheme_names(i) + " scheme");
    
    % part(d): Calculate and plot the error
    [N_list, error, delta_x, delta_t] = calculate_error(a, b, K, T, sigma, q, q_delta, r, g, g1, g2, phi, schemes{1}, i, 'put');
    tabulate(N_list, error, delta_x, delta_t, i);
end

function [h, k, a, b, K, T, sigma, q, q_delta, r, g, g1, g2, phi] = input_pde()
    % This is a helper function to take the input values.
    
    [h, k, a, b, T] = deal(0.01, 0.01, -2, 1, 1);   % a = x_min,    b = x_max
    [r, sigma, delta, K] = deal(0.25, 0.6, 0.2, 10);
    
    q = 2 * r / (sigma ^ 2);
    q_delta = 2 * (r - delta) / (sigma ^ 2);
    
    g_helper = @(x, t, q, q_delta) exp(x * (q_delta - 1) / 2 ) - exp(x * (q_delta + 1) / 2);
    g = @(x, t, q, q_delta) exp(t * ((q_delta - 1) ^ 2 + 4 * q ) / 4) * (g_helper(x, t, q, q_delta) .* (g_helper(x, t, q, q_delta) >= 0));
   
    g1 = @(x, t, q_delta) exp(0.5 * x * (q_delta - 1) + 0.25 * t * (q_delta - 1) ^ 2);
    g2 = @(x, t, q_delta) zeros(size(t));
    
    phi_helper = @(x, q_delta) exp(0.5 * x * (q_delta - 1)) - exp(0.5 * x * (q_delta + 1)); 
    phi = @(x, q_delta) phi_helper(x, q_delta) .* (phi_helper(x, q_delta) >= 0);
end