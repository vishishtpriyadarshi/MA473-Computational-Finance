function compute_delta(V, S, dx, option)
    % This function computes and displays the plot of the delta of the specified 
    % option using the solution matrix.
    
    delta = [];
    [~, len] = size(S);
    for i = 1 : len - 1
        delta = [delta ( V(1, i + 1) - V(1, i) ) / dx];
    end
    
    figure();
    plot(S(1 : end - 1), delta);
    xlabel("Stock price (S) at t = 0");
    ylabel("delta");
    title("Plot of delta for " + option + " Option");
end

