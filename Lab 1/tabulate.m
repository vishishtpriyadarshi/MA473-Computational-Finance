function tabulate(N, error, log_error, idx)
    % This function is used to tabulate the error and the order of convergence
    % for each different N values
    
    schemes = ["FTCS", "BTCS", "Crank-Nicolson"];
    fprintf("********  %s Scheme  ********\n", schemes(idx));
    
    T = table((1 : length(N))', N', error', log_error');
    T.Properties.VariableNames = {'SI No.' 'N' 'Max error (En)' 'Order of Convergence'};
    disp(T);
    
    figure();
    loglog(N, error);
    title("Error ( E_{N} ) vs N: " + schemes(idx) + " scheme");
    xlabel("N");
    ylabel("Error ( E_{N} )");
    
    figure();
    loglog(N, log_error);
    title("Order of Convergence ( log_{2}( E_{N} / E_{2N} ) vs N ) for " + schemes(idx) + " scheme");
    xlabel("N");
    ylabel('$log_{2}( \frac{E_{N}}{ E_{2N}} )$', 'Interpreter', 'latex');
end

