close all, clear all, clc;


%a
%Minsta positiva roten är nästan noll.

%Största rot
%Tar bort tal flera gånger så får man roten:
%7.37
global f;

%b
f = @(x) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - 11 * x .* exp(-x);

%c
function root = find_root_secant(x_0, x_1)
    %Import mathematical functions
    global f;

    x_prev = x_0;
    x = x_1;
    error_prev = 1;

    for i = 1:100
        x_delta = f(x)*((x-x_prev)/(f(x)-f(x_prev)));
        
        %Error and convergence testing
        error = abs(x_delta);
        x_prev = x;
        x = x - x_delta;

        relative_error = error / x;

        K = error / (error_prev^((1+sqrt(5))/2));
        error_prev = error;
        disp(["X:", x, "   K:" , K,]);

        %Exit condition
        if relative_error <= 10^(-8)
            disp(["iterations: ", i]);
            break;
        end
    end
    root = [x, relative_error];
end



%B plotting function--------------------------------------------------------------------------------------------
x = linspace(-1, 8, 500);

%Plot function
subplot(2, 2, 1)
plot(x, f(x));
title("f(x)");
grid on;
ylim([-500, 500]);
hold on;
%Plot zero reference line
plot(x, zeros(length(x)));
hold on;

%C Find roots and plot them-----------------------------------------------------------------------------------------
roots = [find_root_secant(0, 1); find_root_secant(6, 7)];
plot(roots, f(roots), "o");

%Plot first root
subplot(2, 2, 2)
x = linspace(roots(1)-0.1, roots(1)+0.1, 500);
plot(x, f(x));
title(["First root   ", roots(1)]);
grid on;
hold on;
plot(roots(1), f(roots(1)), "o");
hold on;
%Plot zero reference line
plot(x, zeros(length(x)));
hold on;

%Plot second root
subplot(2, 2, 3)
x = linspace(roots(2)-0.1, roots(2)+0.1, 500);
plot(x, f(x));
title(["Second root   ", roots(2)]);
grid on;
hold on;
plot(roots(2), f(roots(2)), "o");
hold on;
%Plot zero reference line
plot(x, zeros(length(x)));
hold on;

disp("The roots are:")
disp(["Root 1: ",roots(1), "   Root 2: ", roots(2)])

%D Definition of convergence--------------------------------------------------------------------------------------------
%Konvergensen är gyllene snittet: lim då (n -> inf); abs(error) / (abs(prev_error)^2) = K


%E
%Bestäm konvergenskonstanten.

%x_root = roots(2);

%x_delta = f(x)*((x-x_prev)/(f(x)-f(x_prev)));
%error1 = abs(x_delta);
%error2 = abs(x_delta = f(x)*((x-x_prev)/(f(x)-f(x_prev))))
%K = error2 / (error1^2);

%K





%F och G. Visa rötter och att konvergens var OK. Konvergens ser man i att K-värdet hålls ungefär samma mellan iterationerna.
disp("Root:                 Relative error:");
disp(roots);