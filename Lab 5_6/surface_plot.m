function surface_plot(x, t, U, idx, K, option)
    % This function creates a surface plot for the corresponding inputs.
    % INPUT: x, y, z (representing the data points)
    % OUTPUT: Surface plot
    
    schemes = ["BTCS", "Crank-Nicolson"];
    figure();
    surf(x, t, U);
    title("Surface plot for δx = δt = 0.01: " + schemes(idx) + " scheme");
    xlabel("x");
    ylabel("t");
    zlabel("u(x, t)");
    colorbar;
    
    % Plot the payoff and price of the option at t = 0
    figure();
    plot(x(1:end-1), U(end, 1:end-1));
    hold on;
    payoff = x(1:end-1);
    for i = 1:length(payoff)
        payoff(i) = payoff(i) - K;
        if strcmp(option, 'put')
            payoff(i) = payoff(i) * (-1);
        end
        if payoff(i) < 0
            payoff(i) = 0;
        end
    end
    
    plot(x(1:end-1), payoff);
    title("Price of the option at t = 0 & Payoff: " + schemes(idx) + " scheme");
    xlabel("S");
    ylabel("U(S, 0)");
    legend("Price of option at t = 0", "max(S - K, 0)");
    hold off;
end

