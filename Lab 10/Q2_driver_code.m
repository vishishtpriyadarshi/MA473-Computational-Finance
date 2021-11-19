[mu, sigma, X0, t0, t1] = deal(-10, 1, 0, 0, 4);
[err_1, err_2, order_of_convg_1, order_of_convg_2] = deal([]);

len = 14;
dt = zeros(len, 1);
dt(1) = 1/8;
for i = 2 : len
    dt(i) = dt(i - 1) / 2;
end

disp("Solving the SDE ...");
for i = 1 : len
    [error_1, error_2] = deal([]);
    for j = 1 : 100
        euler_maruyama_soln_1 = euler_maruyama(mu, sigma, X0, dt(i), t0, t1, "Langevin");
        euler_maruyama_soln_2 = euler_maruyama(mu, sigma, X0, dt(i)/2, t0, t1, "Langevin");
        
        milstein_soln_1 = milstein(mu, sigma, X0, dt(i), t0, t1, "Langevin");
        milstein_soln_2 = milstein(mu, sigma, X0, dt(i)/2, t0, t1, "Langevin");

        error_1(end + 1) = norm(euler_maruyama_soln_1 - euler_maruyama_soln_2(1:2:end));
        error_2(end + 1) = norm(milstein_soln_1 - milstein_soln_2(1:2:end));
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