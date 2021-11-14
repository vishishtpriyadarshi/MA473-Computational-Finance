function surface_plot(x, t, U, idx)
    % This function creates a surface plot for the corresponding inputs.
    % INPUT: x, y, z (representing the data points)
    % OUTPUT: Surface plot
    
    schemes = ["FTCS", "BTCS", "Crank-Nicolson"];
    figure();
    surf(x, t, U);
    
    title("Surface plot for dx = 0.01 & dt = 0.001");
    subtitle(schemes(idx));
    xlabel("R");
    ylabel("t");
    zlabel("H(R, t)");
    colorbar;
end

