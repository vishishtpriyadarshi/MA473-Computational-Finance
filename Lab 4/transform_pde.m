function [U, S, t] = transform_pde(V, xi, tau, T)
    % This function is used to reverse the transformation applied to the
    % Black Scholes equation to get the required solution matrix
    
    q = 0.5;
    S = q * xi ./ (1 - xi);
    
    idx = find(S > 4);
    S = S(1 : idx(1));
    t = T - tau;
    U = (S + q) .* V(:, 1 : idx(1));
    flipud(U);
end

