function [x, t, dx, dt, r] = create_grid(N, M, a, b, T)
    % This function creates a grid/mesh and returns the grid parameters.
    % INPUT: N, M, a, b, T, where
    % N = no of data points on x axis from (a, b)
    % M = no of data points on t axis from (0, T]
    % OUTPUT: x, t, dx, dt, r
    
    x = linspace(a, b, N + 1);
    t = linspace(0, T, M + 1);
    
    dx = 1 / N;
    dt = 1 / M;
    r = dt / dx ^ 2;
end

