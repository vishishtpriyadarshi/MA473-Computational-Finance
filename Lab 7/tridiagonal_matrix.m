function A = tridiagonal_matrix(N, R, r, sigma, h, k, f1, f2, fa, var, sgn)
    % This function creates a tri-diagonal matrix.
    
    A = zeros(N + 1);
    [val_1, val_2, val_3] = deal(f1(R, sigma), f2(R, r), fa(r));

    A(1 : N+2 : end) = var - sgn * 2*val_1*k / h^2 + sgn * val_3*k;
    A(2 : N+2 : end) = sgn * (val_1(2 : N+1)*k / h^2 - val_2(2 : N+1)*k / (2*h));
    A(N+2 : N+2 : end) = sgn * (val_1(1 : N)*k / h^2 + val_2(1 : N)*k / (2*h));

    [A(1,1), A(1, 2), A(1, 3)] = deal(1 - sgn * (3*k)/(2*h), sgn * (4*k)/(2*h), sgn * -k/(2*h));
    [A(N+1, N+1), A(N+1, N)] = deal(1, 0);
end

