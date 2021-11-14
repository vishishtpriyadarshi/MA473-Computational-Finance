function val = linear_basis(x_min, x_max, N, i, x)
    h = (x_max - x_min)/N;
	[x1, x2, x3] = deal(x_min + (i-1)*h, x_min + i*h, x_min + (i+1)*h);
    
    if x < x1 || x > x3
        val = 0;
    else
        if i == N || x < x2
            val = phi(1, x_min, x_max, N, i-1, x, 'function-val');
        else
            val = 1 - phi(1, x_min, x_max, N, i, x, 'function-val');
        end
    end
end