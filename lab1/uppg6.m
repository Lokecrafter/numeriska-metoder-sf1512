close all, clear all, clc;

function ret = own_trapz(integrand, intervals_matrix, start_steps_per_interval, iterations)
    steps = start_steps_per_interval;

    for i = 1:iterations
        evaluated_integral = 0;

        %Integrate over every sub-interval
        size_intervals = size(intervals_matrix);
        for interval = 1:(size_intervals(2))
            start_x = intervals_matrix(1, interval);
            end_x = intervals_matrix(2, interval);
            
            step_size = (end_x - start_x)/steps;
            
            x = start_x + step_size * [0:steps];
            y = integrand(x);

            evaluated_integral = evaluated_integral + (step_size * (sum(y) - (y(1) + y(steps+1)) / 2));
            
        end
        ret(i) = evaluated_integral;
        steps = 2 * steps;
    end
    ret=ret';
end

function ret = richardson_extrapolation(result_one_step_size, result_double_step_size)
    ret = result_one_step_size - ((result_one_step_size - result_double_step_size) / 3);
end


function ret = f(x_values)
    ret = zeros(size(x_values));

    base_function = @(t) (1 - exp(-(t./3).^3)) ./ (4*t.^3);

    for i = 1:length(x_values)
        x = x_values(i);
        if and(-0.1 < x, x < 0.1)
            ret(i) = base_function(-0.1) + (base_function(0.1) - base_function(-0.1)) * (x+0.1);
        else
            ret(i) = base_function(x);
        end
    end
end

a= 0; 
b= 22361; 
trapz_interations = 2;

intervals = [0, 15; 15, 200; 200, b]';

result = own_trapz(@f, intervals, 10, trapz_interations);
err = abs(result(end) - result(end - 1));

xx = linspace(0, 10, 100);
plot(xx, f(xx), "o-");

xlim([-1, 20]);

hold on;
plot(xx, zeros(size(xx)));

disp(result);
disp(["Error: ", err]);
disp("Richardson näst sista: ")
disp(richardson_extrapolation(result(1:(end - 1)), result(2:end)));