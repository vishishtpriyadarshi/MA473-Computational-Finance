function val = simpsons(f, x_min, x_max, h, N, i, j, a, b) 
    val = (h / 6) * (f(x_min, x_max, h, N, i, j, a, 1) + 4*f(x_min, x_max, h, N, i, j, (a + b)/2, 0) + f(x_min, x_max, h, N, i, j, b, 0));
end