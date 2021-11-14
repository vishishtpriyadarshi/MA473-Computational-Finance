function U = BTCS(N, M, R, h, k, T, r, sigma, f1, f2, g, fa)
    % This function solves the given general parabolic PDE using BTCS scheme.
    % OUTPUT: U (representing the solution matrix)
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    B = zeros(N + 1, 1);
    U(1, :) = g(R, T);
    
    for i = 2 : M + 1
        % Prepare tri-diagonal matrix
        A = tridiagonal_matrix(N, R, r, sigma, h, k, f1, f2, fa, 1, 1);
        
		B(2 : N) = U(i-1, 2 : N);
		B(1) = U(i-1, 1);

		U(i, :) = A \ B;
    end
    
	U = flipud(U);
end