clear all, close all, clc
global f;
global y0;
global L;
global number_of_steps;
global tolerance;


f = @(x,y) -(1/6 + (pi.*sin(pi.*x))/(1.6 - cos(pi.*x))).*y;
L = 6;
y0 = 2.5;
number_of_steps = 10;
tolerance = 1e-4;

function ret = do_euler(y_prim, x_span, y0, h, tolerance)
    y_end_prev = 0;
    
    
    for i = 1:20
        xx = x_span(1):h:x_span(2);
        yy = zeros(size(xx));
        
        yy(1) = y0;
        
        for n = 1:length(yy) % Längden av x fast vi vill inte göra en 
            % gång för många så det blir length(x)-1. n är vår loop-variabel.
            yy(n+1) = yy(n) + h * y_prim(xx(n), yy(n)); % Beräknar nästa y-värde med 
            % Eulers metod! So simple! <3 :D
        end
        
        ret = yy;
        
        E_trunk = abs(yy(end) - y_end_prev);
        %disp("Euler   E_trunk: " + E_trunk);
        
        if(E_trunk < tolerance)
            if(i >= 10)
                disp("Iterations euler: " + i);
                break;
            end
        end
        
        y_end_prev = yy(end);
        h = h * 0.5;
    end
end


x_span = [0, L];
h = L / number_of_steps;
result = do_euler(f, x_span, y0, h, tolerance);

xx = linspace(x_span(1), x_span(2), length(result));

plot(xx, result);





function ret = V(xx)
    global f;
    global y0;
    global L;
    global number_of_steps;
    global tolerance;

    ret = zeros(size(xx));

    for i = 1:length(xx)
        if(xx(i) == 0)
            ret(i) = 0;
            disp(xx(i) + "      is zero!!!!!!!!!!!!!!")
            continue;
        end

        x_span = [0, xx(i)];
        h = (i) / number_of_steps;
        result = do_euler(f, x_span, y0, h, tolerance);

        h = xx(i) / (length(result) - 1);

        % Area looked as to be about 3/10 of the plot-square (2.5*6) which gives an area of 4.5. 
        % Only calculating the area of the integral with Trapz's method gives the value 4.5647; we trust this value.
        % integrand = result;
        % volume = h * (sum(integrand) - (integrand(1) + integrand(length(integrand)))) * 0.5;
        integrand = result.*result;
        volume = pi * h * (sum(integrand) - (integrand(1) + integrand(length(integrand)))) * 0.5;
        %disp(['Trapz integral (rotational volume) value: ', num2str(volume)])

        ret(i) = volume;
    end
end


% Integral/volume
% ------------------------------------------------

V([2, 6])

hold on;
xxxxx = 0:0.1:6;
plot(xxxxx, V(xxxxx));




V_65_procent = V(6) * 0.65;
hold on;
plot([0, L], [V_65_procent, V_65_procent])


function root = find_root_secant(f, x_0, x_1, tolerance)
    % Import mathematical functions

    x_prev = x_0;
    x = x_1;
    error_prev = 1;

    disp("     x                K                   Relative_error")

    for i = 1:100
        x_delta = f(x)*((x-x_prev)/(f(x)-f(x_prev)));
        
        %Error and convergence testing
        error = abs(x_delta);
        x_prev = x;
        x = x - x_delta;

        relative_error = error / abs(x);

        K = error / (error_prev^((1+sqrt(5))/2));
        error_prev = error;
        disp([x, K, relative_error]);

        %Exit condition
        if relative_error < tolerance
            if i > 2
                disp("Iterations i: " + i);
                break;
            end
        end
    end
    root = [x, relative_error];
end



error_func = @(x) V(x) - V_65_procent;
x_guess_1 = 2;
x_guess_2 = 2.1;
tolerance = 1e-3;
root = find_root_secant(error_func, x_guess_1, x_guess_2, tolerance);

disp("L 65% of volume: " + root(1) + " +/- " + root(2)*100 + " %")
plot(root(1), V(root(1)), "o");