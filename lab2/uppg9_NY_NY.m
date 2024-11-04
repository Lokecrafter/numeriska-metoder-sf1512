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


x_span = [0, L];
h = L / number_of_steps;





function ret = V(xx, tolerance)
    global f;
    global y0;
    global L;
    global number_of_steps;

    ret = zeros(size(xx));

    for i = 1:length(xx)
        if(xx(i) == 0)
            ret(i) = 0;
            disp(xx(i) + "      is zero!!!!!!!!!!!!!!")
            continue;
        end

        
        h = xx(i) / number_of_steps;
        volume_prev = 0;

        for iteration = 1:30   
            
            %Euler
            xxxx = 0:h:xx(i);
            yy = zeros(size(xxxx));
            yy(1) = y0;
            for n = 1:(length(yy) - 1)
                yy(n+1) = yy(n) + h * f(xxxx(n), yy(n));
            end            
            
            %Trapz
            integrand = yy.*yy;
            volume = pi * h * (sum(integrand) - (integrand(1) + integrand(length(integrand)))) * 0.5;
            ret(i) = volume;


            E_trunk = abs(volume - volume_prev);
            if(E_trunk < tolerance)
                if iteration >= 2
                    disp("Volume iterations: " + iteration + "   E_trunk: " + E_trunk);
                    break
                end
            end

            volume_prev = volume;
            h = h * 0.5;
        end
    end
end


% Integral/volume
% ------------------------------------------------

volume1 = V(6, tolerance);
volume2 = V(6, tolerance * 0.1);
disp("Volume: " + volume2 + "   E_trunk: " + abs(volume2 - volume1));

hold on;
xxxxx = 0:0.2:6;
plot(xxxxx, V(xxxxx, tolerance));




V_65_procent = V(6, tolerance * 0.1) * 0.65;
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



error_func = @(x) V(x, tolerance) - V_65_procent;
x_guess_1 = 2;
x_guess_2 = 2.1;
tolerance = 1e-3;
root1 = find_root_secant(error_func, x_guess_1, x_guess_2, tolerance);
tolerance = tolerance * 0.1;
root2 = find_root_secant(error_func, x_guess_1, x_guess_2, tolerance);

disp("L 65% of volume 1: " + root1(1) + " +/- " + root1(2)*100 + " %")
disp("L 65% of volume 2: " + root2(1) + " +/- " + root2(2)*100 + " %")
disp("L 65% of volume  : " + root2(1) + " +/- " + abs(root1(1) - root2(1)));
plot(root1(1), V(root1(1),tolerance), "o"); hold on;
plot(root2(1), V(root2(1),tolerance), "o");