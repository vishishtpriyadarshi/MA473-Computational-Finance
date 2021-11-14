function [x, tau, h, k] = create_grid(x_min, x_max, N, M, T, sigma)
    x = linspace(x_min, x_max, N + 2);
    tau = linspace(0, T*(sigma^2)/2, M + 1); 
    [h, k] = deal(x(2) - x(1), tau(2) - tau(1));
end