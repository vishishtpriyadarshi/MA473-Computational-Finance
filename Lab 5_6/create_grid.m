function [x, t, lambda, N, M] = create_grid(h, k, a, b, T)
    % This function creates a grid/mesh and returns the grid parameters.
    % INPUT: N, M, a, b, T, where
    % N = no of data points on x axis from (a, b)
    % M = no of data points on t axis from (0, T]
    % OUTPUT: x, t, dx, dt, lambda
    
    N = (b - a) / h;
    M = ceil(T / k);
    
    x = linspace(a, b, N + 1);
    t = linspace(0, T, M + 1);
    
    lambda = k / (h ^ 2);
end

