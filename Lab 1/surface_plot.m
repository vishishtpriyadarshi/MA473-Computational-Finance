function surface_plot(x, t, U, idx)
    % This function creates a surface plot for the corresponding inputs.
    % INPUT: x, y, z (representing the data points)
    % OUTPUT: Surface plot
    
    schemes = ["FTCS", "BTCS", "Crank-Nicolson"];
    figure();
    surf(x, t, U);
    title("Surface plot for N = 10: " + schemes(idx) + " scheme");
    xlabel("x");
    ylabel("t");
    zlabel("u(x, t)");
    colorbar;
end

