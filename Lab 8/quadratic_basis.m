function val = quadratic_basis(x_min, x_max, N, i, x)
    h = (x_max - x_min)/N;
	[x1, x2, x3] = deal(x_min + (i-2)*h, x_min + i*h, x_min + (i+2)*h);
    type = 'function-val';
    
    if x < x1 || x > x3
        val = 0;
    else
        if mod(i, 2) == 1
            val = 4*phi(2, x_min, x_max, N, i-1, x, type) - 4*phi(2, x_min, x_max, N, i-1, x, type)^2;
        elseif i == 0
            if x <= x3
                val = 1 - 3*phi(2, x_min, x_max, N, i, x, type) + 2*phi(2, x_min, x_max, N, i, x, type)^2;
            else
                val = 0;
            end
        elseif i == N
            if x >= x1
                val = -phi(2, x_min, x_max, N, i-2, x, type) + 2*phi(2, x_min, x_max, N, i-2, x, type)^2;
             else
                val = 0;
            end
        else
            if (x1 <= x) && (x <= x2)
                val = - phi(2, x_min, x_max, N, i-2, x, type) + 2*phi(2, x_min, x_max, N, i-2, x, type)^2;
            elseif (x2 <= x) && (x <= x3)
                val = 1 - 3*phi(2, x_min, x_max, N, i, x, type) + 2*phi(2, x_min, x_max, N, i, x, type)^2;
            else
                val = 0;
            end
        end
	end
end

