[mu, sigma, X0, t0, t1] = deal(0.75, 0.3, 307, 0, 1);
[err_1, err_2, order_of_convg_1, order_of_convg_2] = deal([]);

len = 15;
dt = zeros(len, 1);
dt(1) = 1/8;
for i = 2 : len
    dt(i) = dt(i - 1) / 2;
end

disp("Solving the SDE ...");
for i = 1 : len
    [error_1, error_2] = deal([]);
    for j = 1 : 200
        euler_maruyama_soln = euler_maruyama(mu, sigma, X0, dt(i), t0, t1, "Black-Scholes");
        milstein_soln = milstein(mu, sigma, X0, dt(i), t0, t1, "Black-Scholes");
        exact_soln = exact_solution(mu, sigma, X0, dt(i), t0, t1);

        error_1(end + 1) = norm(euler_maruyama_soln - exact_soln);
        error_2(end + 1) = norm(milstein_soln - exact_soln);
    end
    err_1(end + 1) = mean(error_1);
    err_2(end + 1) = mean(error_2);
end

for i = 1 : length(dt) - 1
    order_of_convg_1(end + 1) = log2(err_1(i + 1) / err_1(i));
    order_of_convg_2(end + 1) = log2(err_2(i + 1) / err_2(i));
end

disp("Plotting the numerical solution ...");
figure();
loglog(dt(2:end), order_of_convg_1);
hold on;
loglog(dt(2:end), order_of_convg_2);
hold off;
legend("Euler-Maruyama", "Milstein");
xlabel("Δt");
ylabel("Order of COnvergence");
title("Plot for order of convergence (loglog plot)");