function [N_list, error, log_error] = calculate_error(a, b, T, f, g1, g2, phi, scheme, fact)
    % This function is used to tabulate the error and the order of convergence
    % for each different N values
    % INPUT: a, b, T, f, g1, g2, phi, scheme, idx
    % OUTPUT: N_list, error, log_error
    
    N = 10;
    [N_list, error, log_error] = deal([]);
    
    while N < 160
        [x, t, dx, dt, r] = create_grid(N, 2 * (N^2), a, b, T);
        U1 = scheme(N, 2 * (N^2), x, t, r, dx, dt, f, g1, g2, phi); 
        
        [x, t, dx, dt, r] = create_grid(2*N, 2 * ((2*N) ^ 2), a, b, T);
        U2 = scheme(2*N, 2 * ((2*N) ^ 2), x, t, r, dx, dt, f, g1, g2, phi);
        
        % Pruning the U2 matrix for comparison
        row_step = fix((size(U2, 1) - 1) / (size(U1, 1) - 1));
        col_step = fix((size(U2, 2) - 1) / (size(U1, 2) - 1));
        U2 = U2(1:row_step:end, 1:col_step:end);
        
        error(end + 1) = max(max(abs(U1 - U2)));
        N_list(end + 1) = N;
        N = N * 2;
    end
    
    for i = 1 : length(N_list) - 1
        log_error(end + 1) = (log2(error(i) / error(i + 1))) / (1 / fact);
    end
    
    log_error(end + 1) = log_error(end);
end