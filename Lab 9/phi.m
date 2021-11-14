function val = phi(degree, x_min, x_max, h, N, i, x, type)
	[x1, x2] = deal(x_min + i*h, x_min + (i + degree)*h); 
	
    if x >= x1 && x <= x2
        if strcmp(type, 'derivative')
            val = 1 / h;
        else
            val = (x - x1) / h;
        end
	else
		val = 0;
    end
end