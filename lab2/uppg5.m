% lab 2 upg 5

%-----------------------------------------------

% DELUPPGIFT A

%-----------------------------------------------

close all; clear all; clc;

f =@(k) (3 .* k) ./ 2 - 0.05 - (k .* sin(pi ./ (10 .* k))) ./ (2 .* pi)   -   5 .* (k + 3).^3;

k_start = -4.1;
k_end = -4.0;
xx = linspace(k_start, k_end, 30000);

tolerance = 1e-8;

plot(xx, f(xx), '-')
hold on
plot([k_start, k_end],[0,0])
hold on

%-----------------------------------------------

% DELUPPGIFT D

%-----------------------------------------------

function root = find_root_secant(f, x_0, x_1, tolerance)
    % Import mathematical functions

    x_prev = x_0;
    x = x_1;
    error_prev = 1;

    disp("     x       K       relative_error")

    for i = 1:100
        x_delta = f(x)*((x-x_prev)/(f(x)-f(x_prev)));
        
        %Error and convergence testing
        error = abs(x_delta);
        x_prev = x;
        x = x - x_delta;

        relative_error = error / x;

        K = error / (error_prev^((1+sqrt(5))/2));
        error_prev = error;
        disp([x, K, relative_error]);

        %Exit condition
        if relative_error < tolerance
            if i > 2
                disp(["Iterations i: " + i]);
                break;
            end
        end
    end
    root = [x, relative_error];
end

root = find_root_secant(f, k_start, k_end, tolerance);

disp("Root k = " + root(1))
disp("Relative error = " + root(2))

plot(root(1), f(root(1)), 'o')