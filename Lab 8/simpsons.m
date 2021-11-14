function val = simpsons(f, x_min, x_max, N, i, j, a, b) 
    h = b - a;
    val = (h / 6) * (f(x_min, x_max, N, i, j, a, 1) + 4*f(x_min, x_max, N, i, j, (a + b)/2, 0) + f(x_min, x_max, N, i, j, b, 0));
end