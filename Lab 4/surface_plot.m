function surface_plot(x, t, U, idx, m)
    % This function creates a surface plot for the corresponding inputs.
    % INPUT: x, y, z (representing the data points)
    % OUTPUT: Surface plot
    
    if m ~= ""
        m = " using " + m;
    end
    
    schemes = ["FTCS", "BTCS", "Crank-Nicolson"];
    figure();
    surf(x, t, U);
    zlim([0 inf])
    title("Surface plot for dx = 0.01");
    subtitle(schemes(idx) + " scheme" + m);
    xlabel("x");
    ylabel("t");
    zlabel("u(x, t)");
    colorbar;
end

