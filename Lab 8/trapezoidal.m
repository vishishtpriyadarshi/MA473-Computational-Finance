function val = trapezoidal(f, x_min, x_max, N, i, j, a, b) 
    h = b - a;
    val = (h / 2) * (f(x_min, x_max, N, i, j, a, 1) + f(x_min, x_max, N, i, j, b, 0));
end