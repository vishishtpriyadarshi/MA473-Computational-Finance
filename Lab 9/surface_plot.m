function surface_plot(x, tau, V, T, sigma, K, idx, plot_title)
    figure();
    S = K * exp(x);
    time = T - tau * 2 / sigma^2;
    
    surf(S', time', V');
    
    xlabel("S");
    ylabel("t");
    zlabel("V(S, t)");
    
    if idx == 1
        type = " Trapezoidal Rule";
    else
        type = " Simpson's Rule";
    end
    
    title(plot_title + type);
end

