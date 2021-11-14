function val = derivative_linear_basis(x_min, x_max, N, i, x, flag)
     h = (x_max - x_min)/N;
	[x1, x2, x3] = deal(x_min + (i-1)*h, x_min + i*h, x_min + (i+1)*h);

    if flag == 0
		x = x - h/3;
	else
		x = x + h/3;
    end
    
    if x < x1 || x > x3
        val = 0;
    else
        if i == N || x < x2
            val = phi(1, x_min, x_max, N, i-1, x, 'derivative');
        else
            val = -phi(1, x_min, x_max, N, i, x, 'derivative');
        end
    end
end

