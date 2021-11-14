function compute_gamma(V, S, dx, option)
    % This function computes and displays the plot of the gamma of the specified 
    % option using the solution matrix.
    
    gamma = [];
    [~, len] = size(S);
    for i = 2 : len - 1
        gamma = [gamma ( V(1, i + 1) - 2 * V(1, i) + V(1, i - 1)) / dx ^ 2];
    end
    
    figure();
    plot(S(2 : end - 1), gamma);
    xlabel("Stock price (S) at t = 0");
    ylabel("gamma");
    title("Plot of gamma for " + option + " Option");
end

