close all, clear all, clc;


function root = find_root_secant(f, x_0, x_1, tolerance)
    %Import mathematical functions

    x_prev = x_0;
    x = x_1;
    error = 0;

    for i = 1:100
        x_delta = f(x)*((x-x_prev)/(f(x)-f(x_prev)));
        
        %Error and convergence testing
        error = abs(x_delta);
        x_prev = x;
        x = x - x_delta;

        %Exit condition
        if error <= tolerance
            break;
        end
    end
    root = [x, error];
end



f = @(x) 197 * exp(-((19*x - pi) / 0.003).^2);


xx = linspace(0.164, 0.166, 10000);

plot(xx, f(xx));
hold on;
grid on;
plot([0.164, 0.166], [0.01, 0.01]);
ylim([0, 0.2]);

%Kollar var f = 0.01 grafiskt
%Får:

x_start_guess = [0.1648, 0.1649];
x_end_guess = [0.1658, 0.1659];

%close all;
hold on;

secant_expression = @(x) f(x) - 0.01;

x_start = find_root_secant(secant_expression, x_start_guess(1), x_start_guess(2), 10^(-10));
x_end = find_root_secant(secant_expression, x_end_guess(1), x_end_guess(2), 10^(-10));

disp(["Start:", x_start,"End:", x_end]);

plot(x_start, f(x_start), "o");
hold on;
plot(x_end, f(x_end), "o");

a = x_start(1);
b = x_end(1);
individual_tolerance = 10^(-10);

result1 = quad(f, 0, a, individual_tolerance);
result2 = quad(f, a, b, individual_tolerance);
result3 = quad(f, b, b + a, individual_tolerance);
result4 = quad(f, b, 6, individual_tolerance);

result_quad = result4 + result1 + result3 + result2;
err_quad = 4 * individual_tolerance;

disp(["Result quad:", result_quad, "Error (+/-):", err_quad]);

result5 = integral(f, 0, a, 'AbsTol', individual_tolerance);
result6 = integral(f, a, b, 'AbsTol', individual_tolerance);
result7 = integral(f, b, b + a, 'AbsTol', individual_tolerance);
result8 = integral(f, b, 6, 'AbsTol', individual_tolerance);

result_integral = result8 + result7 + result6 + result5;
err_integral = 4 * individual_tolerance;

disp(["Result integral:", result_integral, "Error (+/-):", err_integral]);

