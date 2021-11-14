function [V, x, t] = transform_pde(y, x, tau, K, T, sigma, q, q_delta)
    % This function is used to reverse the transformation applied to the
    % Black Scholes equation to get the required solution matrix
    
    V = zeros(size(y));
    for i = 1 : length(tau)
        for j = 1 : length(x)
            V(i, j) = K * exp(-0.5 * (q_delta - 1) * x(j) - (0.25 * ((q_delta - 1) ^ 2) + q) * tau(i)) * y(i, j);
        end
    end

    x = K * exp(x);
    t = T - 2 * tau / (sigma ^ 2);
end

