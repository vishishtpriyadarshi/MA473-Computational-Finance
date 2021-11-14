function [x, t, N, M] = create_grid(dx, dt, x_min, x_max, T)
    % This function creates a grid/mesh and returns the grid parameters.
    % INPUT: dx, dt, x_min, x_max, T
    % OUTPUT: x, t, N, M
    % where,
    % N = no of data points on x axis from (x_min, x_max)
    % M = no of data points on t axis from (0, T]
    
    N = ceil((x_max - x_min) / dx);
    M = ceil(T / dt);
    
    x = x_min: dx: x_max;
    t = 0: dt: T;
end

