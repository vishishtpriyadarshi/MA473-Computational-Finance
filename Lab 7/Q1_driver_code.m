schemes = { @(N, M, x, dx, dt, T, r, sigma, f1, f2, g, fa) FTCS(N, M, x, dx, dt, T, r, sigma, f1, f2, g, fa),
            @(N, M, x, dx, dt, T, r, sigma, f1, f2, g, fa) BTCS(N, M, x, dx, dt, T, r, sigma, f1, f2, g, fa),
            @(N, M, x, dx, dt, T, r, sigma, f1, f2, g, fa) CrankNicolson(N, M, x, dx, dt, T, r, sigma, f1, f2, g, fa)};
      
for i = 1 : length(schemes)
    % Take input variables for the PDE
    [dx, dt, x_min, x_max, T, K, r, sigma, f1, f2, g, fa] = input_pde();
    
    % Create grid and get mesh parameters
    [x, t, N, M] = create_grid(dx, dt, x_min, x_max, T);
    
    % Apply suitable scheme to get the solution matrix
    V = schemes{i}(N, M, x, dx, -dt, T, r, sigma, f1, f2, g, fa);

    % Display the surface plot
    surface_plot(x, t, V, i);
end


function [dx, dt, x_min, x_max, T, K, r, sigma, f1, f2, g, fa] = input_pde()
    % This is a helper function to take the input values.
    
    dx = 0.01;
    dt = 0.001;
     
    [x_min, x_max, T, K] = deal(0, 1, 0.2, 100);
    [r, sigma] = deal(0.05, 0.25);
    
    f1 = @(R, sigma) 0.5 * (sigma * R) .^ 2;
    f2 = @(R, r) 1 - r * R;  
    
    helper = @(R, T) 1 - R / T; 
    g = @(R, T) helper(R, T) .* (helper(R, T) >= 0);
    
    fa = @(R) zeros(size(R));
end