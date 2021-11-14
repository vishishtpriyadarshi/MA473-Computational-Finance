function val = phi(degree, x_min, x_max, N, i, x, type)
	h = (x_max - x_min)/N;
	[x1, x2] = deal(x_min + i*h, x_min + (i + degree)*h); 
	
    if x >= x1 && x <= x2
        if strcmp(type, 'derivative')
            val = 1 / h;
        else
            val = (x - x1) / (x2 - x1);
        end
	else
		val = 0;
    end
end