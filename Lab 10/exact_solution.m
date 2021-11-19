function X = exact_solution(mu, sigma, X0, dt, t0, t1)
    t = t0 : dt : t1;
    X = zeros(size(t));
    n = length(t);
    
    X(1) = X0;
    W = sqrt(dt) * randn(n, 1);
    
    for i = 2 : n
        X(i) = X(i - 1) * exp((mu - 0.5 * sigma ^ 2) * dt + sigma * W(i));
    end
end