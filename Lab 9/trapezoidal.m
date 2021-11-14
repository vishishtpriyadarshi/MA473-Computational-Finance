function val = trapezoidal(f, x_min, x_max, h, N, i, j, a, b) 
    val = (h / 2) * (f(x_min, x_max, h, N, i, j, a, 1) + f(x_min, x_max, h, N, i, j, b, 0));
end