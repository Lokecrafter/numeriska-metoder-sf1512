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
        if and(-0.01 < x, x < 0.01)
            ret(i) = base_function(-0.01) + (base_function(0.01) - base_function(-0.01)) * (x+0.01);
        else
            ret(i) = base_function(x);
        end
    end
end

a= 0; 
b= 35356; 
trapz_per_interval = 10;
trapz_iterations = 12;

intervals = [5000, b; 200, 5000; 15, 200; a, 15]'; %Ombytt för 

richard = own_trapz(@f, intervals, trapz_per_interval, trapz_iterations);

tail_cut_err = 10^(-10);
begin_approx_err = 1.1 * 10^(-11);
err = begin_approx_err + tail_cut_err + abs(richard(2:end) - richard(1:(end - 1)));
%Felet från approximationen i det "svajiga" området för gränsvärdet då x går mot 0 blir exakt 1.1e-11 vilket säger att vår approximation är tillräckligt bra för att inte påverka resultatet.


disp(richard);
disp("Error: ");
disp(err);
richard = richardson_extrapolation(richard(1:(end - 1)), richard(2:end));
disp("Richardson extrapolation: ")
disp(richard);
disp(["Result: ", richard(end), " +/- ", err(end)]);




xx = linspace(0, 10, 100);
plot(xx, f(xx), "o-");

xlim([-1, 20]);

hold on;
plot(xx, zeros(size(xx)));