function A = tridiagonal_matrix(N, a, b, c)
    % This function creates a NxN tri-diagonal matrix.
    % INPUT: N, a, b, c, where a = values on the main diagonal, b = values
    % on the first diagonal above the main and c = values on the diagonal
    % below the main diagonal.
    % OUTPUT: Tri-diagonal matrix
    
    A = diag(a * ones(1, N)) + diag(b * ones(1, N - 1), 1) + diag(c * ones(1, N - 1), -1);
end

