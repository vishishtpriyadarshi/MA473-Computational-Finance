function U = CrankNicolson(N, M, R, h, k, T, r, sigma, f1, f2, g, fa)
    % This function solves the given general parabolic PDE using CrankNicolson scheme.
    % OUTPUT: U (representing the solution matrix)
    
    % Create solution matrix with given ICs
    U = zeros(M + 1, N + 1);
    B = zeros(N + 1, 1);
    U(1, :) = g(R, T);
    
    for i = 2 : M + 1
        % Prepare tri-diagonal matrix
        A = tridiagonal_matrix(N, R, r, sigma, h, k, f1, f2, fa, 2, 1);
        
        [val_1, val_2, val_3] = deal(f1(R, sigma), f2(R, r), fa(r));
		B(2 : N) = (-val_1(2 : N)*k/h^2 + val_2(2 : N)*k/(2*h)) .* U(i-1,1 : N - 1) + (2 + 2*val_1(2 : N)*k/h^2 - val_3*k) .* U(i-1, 2 : N) ...
                    + (-val_1(2 : N)*k/h^2 - val_2(2 : N)*k/(2*h)) .* U(i-1,3 : N+1);
		B(1) = U(i-1, 1);

		U(i, :) = A \ B;
    end
    
	U = flipud(U);
end


