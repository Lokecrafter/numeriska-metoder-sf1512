close all, clear all, clc;




function ret = integrand(x)
    base_function = @(x) (1 - exp(-(x/3).^3)) ./ (4*x.^3);
    if and((-0.1 < x), (x < 0.1))
        ret = base_function(-0.1) + (base_function(0.1) - base_function(-0.1)) * (x+0.1);
    else
        ret = base_function(x);
    end
end

a= 0; 
middle_point = 100;
b= 22361; 
divide_step_size_iterations = 10;

interval1 = [a, 0.1];

interval2 = [0.1, middle_point];

interval3 = [middle_point, b];


result = 0;

result = result + own_trapz(integrand3,interval3,10,divide_step_size_iterations);
result = result + own_trapz(integrand2,interval2,2,divide_step_size_iterations); %Short interval. Not many splits needed
result = result + own_trapz(integrand1,interval1,1,1); %Linear interpolation. Only one trapz needed.

err = abs(result(end) - result(end - 1));

xx = [linspace(interval1(1),interval1(2)), linspace(interval2(1),interval2(2)), linspace(interval3(1),interval3(2))];
plot(xx, integrand(xx), "o-");


disp(result);
disp(["Error: ", err]);
disp("Richardson näst sista: ")
disp(richardson_extrapolation(result(1:(end - 1)), result(2:end), 2));