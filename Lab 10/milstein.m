function X = milstein(mu, sigma, X0, dt, t0, t1, type)
    t = t0 : dt : t1;
    X = zeros(size(t));
    n = length(t);
    
    X(1) = X0;
    W = sqrt(dt) * randn(n, 1);
    
    for i = 2 : n
        if type == "Langevin"
            X(i) = X(i - 1) + mu * X(i - 1) * dt + sigma * W(i);
        else
            X(i) = X(i - 1) + mu * X(i - 1) * dt + sigma * X(i - 1) * W(i) + ...
                    0.5 * (sigma * X(i - 1)) * sigma * (W(i) ^ 2 - dt);
        end
    end
end