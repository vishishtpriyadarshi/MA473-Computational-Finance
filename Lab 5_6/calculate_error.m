function [N_list, error, delta_x, delta_t] = calculate_error(a, b, K, T, sigma, q, q_delta, r, g, g1, g2, phi, scheme, i, option)
    % This function is used to tabulate the error and the order of convergence
    % for each different delta_x (or N) values
    
    if i == 1
        [N, M] = deal(4);
    else
        [N, M] = deal(5);
    end
    [N_list, error, delta_x, delta_t] = deal([]);
    fl = 0;
    
    while N <= 320
        fprintf("Calculating error for N = %d\n", N);
        % To speed-up the processing
        if fl == 0
            [x, tau, lambda, h, k] = create_grid_modified(N, M, a, b, T * (sigma ^ 2) / 2);
            [U1, ~, ~] = scheme(N, M, x, tau, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, option);
        end
  
        [x, tau, lambda, h, k] = create_grid_modified(2*N, 2*M, a, b, T * (sigma ^ 2) / 2);
        [U2, ~, ~] = scheme(2*N, 2*M, x, tau, lambda, K, T, sigma, q, q_delta, g, g1, g2, phi, option);
        
        % Pruning the U2 matrix for comparison
        U2_pruned = U2(1:2:end, 1:2:end);
        
        error(end + 1) = max(max(abs(U1 - U2_pruned)));
        N_list(end + 1) = N;
        [delta_x(end + 1), delta_t(end + 1)] = deal(2*h, 2*k);
        [N, M] = deal(N * 2);
        
        % To speed-up the processing
        [U1, fl] = deal(U2, 1);
    end
end

function [x, t, lambda, h, k] = create_grid_modified(N, M, a, b, T)
    h = (b - a) / N;
    k = T / M;
    
    x = linspace(a, b, N + 1);
    t = linspace(0, T, M + 1);
    
    lambda = k / (h ^ 2);
end