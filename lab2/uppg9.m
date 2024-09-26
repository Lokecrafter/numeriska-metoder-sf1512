<<<<<<< Updated upstream
close all, clear all, clc;

=======
% Uppgift 9 - Numerisk integration: Rotationssymmetrisk lur


% Defining function f = y prim
>>>>>>> Stashed changes
f = @(y,x) -(1/6 + (pi.*sin(pi.*x))/(1.6 - cos(pi.*x))).*y;



% Function for Euler's method
% ------------------------------------------------

function ret = do_euler(f, start_x, start_y, step_size, end_x, tolerance)

    prev_end_y = 0;

    for i = 1:100
        x = start_x:step_size:end_x;
        y = zeros(size(x));
        y(1) = start_y;

        % disp(["Number of steps: ", length(x)]);

        for n = 1:length(y)
            y(n+1) = y(n) + step_size * f(y(n), x(n));
        end

        ret = y;
        
        
        error = abs(prev_end_y - y(end));

        % Exit condition based on tolerance.
        if error < tolerance
            if(i ~= 1)
                break;
            end
        end
        
        prev_end_y = y(end);
        step_size = step_size * 0.5; % 0.5 = 1/2. Datorn räknar multiplikation snabbare än division
    end
end


% Run program
% ------------------------------------------------

result = do_euler(f, 0, 2.5, 0.5, 6, 10^(-4));
x = linspace(0,6,length(result));
step_size = 6/(length(result)-1); % Interval of x [0,6] divided by the length of 'result's (i.e number of steps)
disp(['Last y-term: ', num2str(result(end))])
disp(['Step size: ', num2str(step_size)])


% Plot of do_euler
% ------------------------------------------------

plot(x,result,'-o'); 
% '-o' : 
% '-' Skapar en solid linje som förbinder datapunkterna. 
% 'o': Lägger till cirkelmarkörer vid varje datapunkt.
xlabel('x');
ylabel('y');
grid on;
title('Eulers metod')

% Integral
% ------------------------------------------------

integrand = result.*result;
volume = pi * step_size * (sum(integrand) - (integrand(1) + integrand(length(integrand))) * 0.5);
disp(['Trapz integral value: ', num2str(volume)])