quadrature_rules = { @(f, x_min, x_max, h, N, i, j, a, b) trapezoidal(f, x_min, x_max, h, N, i, j, a, b),
                     @(f, x_min, x_max, h, N, i, j, a, b) simpsons(f, x_min, x_max, h, N, i, j, a, b) };
                 
basis_functions = { @(x_min, x_max, h, N, i, x) linear_basis(x_min, x_max, h, N, i, x) };

basis_derivatives = { @(x_min, x_max, h, N, i, x, flag) derivative_linear_basis(x_min, x_max, h, N, i, x, flag) };


[T, K] = deal(1, 10);
[r, sigma, delta] = deal(0.06, 0.3, 0);

q = 2 * r / (sigma ^ 2);
q_delta = 2 * (r - delta) / (sigma ^ 2);

[x_min, x_max, N, M] = deal(-4, 2, 80, 50);
[x, tau, h, k] = create_grid(x_min, x_max, N, M, T, sigma);
x_eff = x(2 : length(x)-1); 

for q_idx = 1 : length(quadrature_rules)
    if q_idx == 1
        fprintf("(a) Solving European Put option using linear basis & Trapezoidal Rule ...\n\n");
    else
        fprintf("(b) Solving European Put option using linear basis & Simpson's Rule ...\n\n");
    end
    
    [quadrature, basis, der_b] = deal(quadrature_rules{q_idx}, basis_functions{1}, basis_derivatives{1});

    [A, B] = deal(zeros(length(x_eff)));
    for i = 1 : length(x_eff)
        for j = 1 : length(x_eff)
            if i == j
                for idx = 0 : 1
                    func1 = @(x_min, x_max, h, N, i, j, val, fl) der_b(x_min, x_max, h, N, i, val, fl) * der_b(x_min, x_max, h, N, j, val, fl);
                    A(i, j) = A(i, j) + quadrature(func1, x_min, x_max, h, length(x), i, j, x(i+idx), x(i+idx+1));

                    func3 = @(x_min, x_max, h, N, i, j, val, fl) basis(x_min, x_max, h, N, i, val) * basis(x_min, x_max, h, N, j, val);
                    B(i, j) = B(i, j) + quadrature(func3, x_min, x_max, h, length(x), i, j, x(i+idx), x(i+idx+1));
                end
            elseif abs(i - j) == 1
                idx = min(i, j);
                func1 = @(x_min, x_max, h, N, i, j, val, fl) der_b(x_min, x_max, h, N, i, val, fl) * der_b(x_min, x_max, h, N, j, val, fl);
                A(i, j) = A(i, j) + quadrature(func1, x_min, x_max, h, length(x), i, j, x(idx+1), x(idx+2));

                func3 = @(x_min, x_max, h, N, i, j, val, fl) basis(x_min, x_max, h, N, i, val) * basis(x_min, x_max, h, N, j, val);
                B(i, j) = B(i, j) + quadrature(func3, x_min, x_max, h, length(x), i, j, x(idx+1), x(idx+2));
            end
        end
    end
    B_final = B + (0.5 * k) .* A;
    A_final = B - (0.5 * k) .* A;

    [w, b] = deal(zeros(length(x_eff), length(tau)));

    for i = 1 : length(tau)
        b(:, i) = ((x_eff - x_min)/(x_max - x_min)) * h * (0.25 * (q_delta - 1) ^ 2) * alpha(x_max, tau(i), q_delta) + (0.25 * (q_delta - 1) ^ 2) * alpha(x_max, tau(i), q_delta);
    end

    % Construct weights
    w(: , 1) = IC(x_eff, q_delta) - phi_b(x_eff, 0, q_delta, x_min, x_max);
    for i = 2 : length(tau)
        w(:, i) = B_final \ (A_final * w(:, i - 1) - (k/2) * (b(:, i) + b(:, i-1)));
    end

    % Get option price at nodes
    base = zeros(N + 1, M + 1);
    for i = 1 : length(x)
        for j = 1 : length(tau)
            base(i, j) = phi_b(x(i), tau(j), q_delta, x_min, x_max);
        end
    end

    temp = zeros(1, M + 1);
    nodes = [temp; w; temp];
    nodes = nodes + base;

    V = nodes;
    temp1 = (-0.5) * (q - 1) * x;
    temp2 = (-0.25) * ((q + 1)^2) * tau;

    for i = 1 : length(x)
        for j = 1 : length(tau)
            V(i, j)=(K * exp(temp1(i) + temp2(j))) * nodes(i, j);
        end
    end

    % Plot graphs
    surface_plot(x, tau, V, T, sigma, K, q_idx, "European Put option using");
end


% ***************  HELPER FUNCTIONS  ***************

function y = IC(x, q_delta)
    y = max(0, exp(0.5 * x * (q_delta - 1)) - exp(0.5 * x * (q_delta + 1)));
end

% Boundary condition at x_min
function y = alpha(x_min, tau, q_delta)
    y = exp(0.5 * (q_delta - 1) * x_min + 0.25 * tau * (q_delta - 1) ^ 2);
end

% Boundary condition at x_max
function y = beta(x_max, tau, q_delta)
    y = zeros(size(tau));
end

function val = phi_b(x, tau, q_delta, x_min, x_max)
    val = beta(x_max, tau, q_delta) - alpha(x_min, tau, q_delta);
    val = alpha(x_min, tau, q_delta) + val * (x - x_min) / (x_max - x_min);
end