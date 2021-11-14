schemes = { @(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g1, g2, phi) FTCS(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g1, g2, phi);
            @(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g1, g2, phi) BTCS(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g1, g2, phi);
            @(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g1, g2, phi) CrankNicolson(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g1, g2, phi)}; 

for i = 1 : length(schemes)
    % Take input variables for the PDE
    [h, k, a, b, K, T, sigma, q, q_delta, r, g1, g2, phi] = input_pde();
    
    % Create grid and get mesh parameters
    [x, tau, lambda, N, M] = create_grid(h, k, a, b, T * (sigma ^ 2) / 2);
    
    % Apply suitable scheme to get the solution matrix
    [V, S, t] = schemes{i}(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g1, g2, phi);
    
    % Display the surface plot
    surface_plot(S, t, V, i);
end

function [h, k, a, b, K, T, sigma, q, q_delta, r, g1, g2, phi] = input_pde()
    % This is a helper function to take the input values.
    
    h = 0.05;
    k = h^2 / 2;
    a = -4;
    b = 2;
    T = 1;
    
    r = 0.06;
    sigma = 0.3;
    delta = 0;
    K = 10;
    
    q = 2 * r / (sigma ^ 2);
    q_delta = 2 * (r - delta) / (sigma ^ 2);
    
    g1 = @(x, t) exp(0.5 * x * (q_delta - 1) + 0.25 * t * (q_delta - 1) ^ 2);
    g2 = @(x, t) zeros(size(t));
    
    phi_helper = @(x) exp(0.5 * x * (q_delta - 1)) - exp(0.5 * x * (q_delta + 1)); 
    phi = @(x) phi_helper(x) .* (phi_helper(x) >= 0);
end