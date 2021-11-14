schemes = { @(N, M, x, t, lambda, dx, dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m) FTCS(N, M, x, t, lambda, dx, dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m),
            @(N, M, x, t, lambda, dx, dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m) BTCS(N, M, x, t, lambda, dx, dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m),
            @(N, M, x, t, lambda, dx, dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m) CrankNicolson(N, M, x, t, lambda, dx, dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m)};

for i = 1 : length(schemes)
    % Take input variables for the PDE
    [dx, dt, x_min, x_max, T, K, r, sigma, delta, f1, f2, g, a, b, c, d] = input_pde();
    
    % Create grid and get mesh parameters
    [x, t, N, M, lambda] = create_grid(dx, dt, x_min, x_max, T);
    
    if i == 1
        methods = [""];
    else
        methods = ["right-matrix division", "Jacobi method", "Gauss-Seidel method", "SOR method"];
    end
    
    for m = methods
        % Apply suitable scheme to get the solution matrix
        U = schemes{i}(N, M, x, flip(t), lambda, dx, -dt, K, T, r, sigma, delta, f1, f2, g, a, b, c, d, m);

        % Display the surface plot
        surface_plot(x, t, U, i, m);
    end
end

function [dx, dt, x_min, x_max, T, K, r, sigma, delta, f1, f2, g, a, b, c, d] = input_pde()
    % This is a helper function to take the input values.
    
    dx = 0.75;
    dt = dx^2 / 80;
    
    [x_min, x_max, T] = deal(0, 30, 1);
    [r, sigma, delta, K] = deal(0.06, 0.3, 0, 10);
    
    f1 = @(x, t, K, r, T) K * exp(-r * (T - t)) - x;  
    f2 = @(x, t, K, r, T) zeros(size(t));
    
    helper = @(x, K) K - x; 
    g = @(x, K) helper(x, K) .* (helper(x, K) >= 0);
    
    a = @(x, t) 1;
    b = @(x, t, sigma) 0.5 * (sigma ^ 2) * (x .^ 2);
    c = @(x, t, r, delta) (r - delta) .* x;
    d = @(x, t, r) -r;
end