schemes = { @(N, M, x, t, lambda, dx, dt, T, r, sigma, delta, f1, f2, g, a, b, c, d, m) FTCS(N, M, x, t, lambda, dx, dt, T, r, sigma, delta, f1, f2, g, a, b, c, d, m),
            @(N, M, x, t, lambda, dx, dt, T, r, sigma, delta, f1, f2, g, a, b, c, d, m) BTCS(N, M, x, t, lambda, dx, dt, T, r, sigma, delta, f1, f2, g, a, b, c, d, m),
            @(N, M, x, t, lambda, dx, dt, T, r, sigma, delta, f1, f2, g, a, b, c, d, m) CrankNicolson(N, M, x, t, lambda, dx, dt, T, r, sigma, delta, f1, f2, g, a, b, c, d, m)};

for i = 1 : length(schemes)
    % Take input variables for the PDE
    [dx, dt, x_min, x_max, T, r, sigma, delta, f1, f2, g, a, b, c, d] = input_pde();
    
    % Create grid and get mesh parameters
    [x, t, N, M, lambda] = create_grid(dx, dt, x_min, x_max, T);
    
    if i == 1
        methods = [""];
    else
        methods = ["right-matrix division", "Jacobi method", "Gauss-Seidel method", "SOR method", "Conjugate Gradient method"];
    end
    
    for m = methods
        % Apply suitable scheme to get the solution matrix
        [V, V_untranformed, S, S_untransformed, time] = schemes{i}(N, M, x, t, lambda*dx, dx, dt, T, r, sigma, delta, f1, f2, g, a, b, c, d, m);

        % Display the surface plot
        surface_plot(S, time, V, i, m);
    end
end

% Compute Greeks
compute_delta(V_untranformed, S_untransformed, dx, "Put");
compute_gamma(V_untranformed, S_untransformed, dx, "Put");


function [dx, dt, x_min, x_max, T, r, sigma, delta, f1, f2, g, a, b, c, d] = input_pde()
    % This is a helper function to take the input values.
    
    dx = 0.01;
    dt = dx / 2;
     
    [x_min, x_max, T] = deal(0, 1, 1);
    [r, sigma, delta] = deal(0.04, 0.25, 0.1);
    
    f1 = @(x, t, r, T) exp(-delta * t);  
    f2 = @(x, t, r, T) zeros(size(t));
    
    helper = @(x) 1 - 2 * x; 
    g = @(x) helper(x) .* (helper(x) >= 0);
    
    a = @(x, t) 1;
    b = @(x, t, sigma) -0.5 * (sigma ^ 2) * (x .^ 2) * (1 - x .^ 2);
    c = @(x, t, r, delta) -(r - delta) .* x .* (1 - x);
    d = @(x, t, r, delta) r * (1 - x) + delta * x;
end