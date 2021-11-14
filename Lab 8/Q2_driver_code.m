quadrature_rules = { @(f, x_min, x_max, N, i, j, a, b) trapezoidal(f, x_min, x_max, N, i, j, a, b),
                     @(f, x_min, x_max, N, i, j, a, b) simpsons(f, x_min, x_max, N, i, j, a, b) };
                 
basis_functions = { @(x_min, x_max, N, i, x) linear_basis(x_min, x_max, N, i, x),
                    @(x_min, x_max, N, i, x) quadratic_basis(x_min, x_max, N, i, x) };

basis_derivatives = { @(x_min, x_max, N, i, x, flag) derivative_linear_basis(x_min, x_max, N, i, x, flag), 
                      @(x_min, x_max, N, i, x, flag) derivative_quadratic_basis(x_min, x_max, N, i, x, flag) };
          
                  
for b_idx = 1 : length(basis_functions)
    for q_idx = 1 : length(quadrature_rules)
        fprintf("\nSolving BVP using %s\n\n", disp_helper(b_idx, q_idx));
        
        [basis, quadrature, der_b] = get_required_functions(quadrature_rules, basis_functions, basis_derivatives, b_idx, q_idx);
        [x_min, x_max, N] = deal(0, 1, 100);
        h = (x_max - x_min)/N;
        
        X = linspace(x_min, x_max, N+1);
        [A_1, A_2, A_3] = deal(zeros(N+1, N+1));
        d = zeros(N+1, 1);
        
        for i = 1 : N - 1
            for j = 0 : N
                for k = 0 : N - 1
                    func1 = @(x_min, x_max, N, i, j, val, fl) der_b(x_min, x_max, N, i, val, fl) * der_b(x_min, x_max, N, j, val, fl);
                    A_1(i+1, j+1) = A_1(i+1, j+1) + quadrature(func1, x_min, x_max, N, i, j, X(k+1), X(k+2));
                    
                    func2 = @(x_min, x_max, N, i, j, val, fl) p(val) * basis(x_min, x_max, N, i, val) * der_b(x_min, x_max, N, j, val, fl);
                    A_2(i+1, j+1) = A_2(i+1, j+1) + quadrature(func2, x_min, x_max, N, i, j, X(k+1), X(k+2));

                    func3 = @(x_min, x_max, N, i, j, val, fl) q(val) * basis(x_min, x_max, N, i, val) * basis(x_min, x_max, N, j, val);
                    A_3(i+1, j+1) = A_3(i+1, j+1) + quadrature(func3, x_min, x_max, N, i, j, X(k+1), X(k+2));

                    if j == 0
                        func4 = @(x_min, x_max, N, i, j, val, fl) f(val) * basis(x_min, x_max, N, i, val);
                        d(i + 1) = d(i + 1) + quadrature(func4, x_min, x_max, N, i, j, X(k+1), X(k+2));
                    end
                end
            end
        end
        
        A = A_1 + A_2 + A_3;
        [A(1, 1), A(N + 1, N + 1)] = deal(1);
        U = A \ d;
        
        % Plot the solution
        figure();
        plot(X, U');
        xlabel("x");
        ylabel("u");
        title(disp_helper(b_idx, q_idx));
    end
end


function val = p(x)
    val = zeros(size(x));
end

function val = q(x)
    val = 2*x + 1;
end

function val = f(x)
    val = sin(x);
end

function [basis, quadrature, der_b] = get_required_functions(quadrature_rules, basis_functions, basis_derivatives, b_idx, q_idx)
    [quadrature, basis, der_b] = deal(quadrature_rules{q_idx}, basis_functions{xor(b_idx, b_idx) + 1}, basis_derivatives{xor(b_idx, b_idx) + 1});
end

function str = disp_helper(i, j)
    if i == 1 && j == 1
        str = "Linear basis function and Trapezoidal rule";
    elseif i == 1 && j == 2
        str = "Linear basis function and Simpson's rule";
    elseif i == 2 && j == 1
        str = "Quadratic basis function and Trapezoidal rule";
    else
        str = "Quadratic basis function and Simpson's rule";
    end
end