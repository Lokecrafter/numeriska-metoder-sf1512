close all, clear all, clc;


%a
%Minsta positiva roten är nästan noll.

%Största rot
%Tar bort tal flera gånger så får man roten:
%7.37

global f;
global f_prim;
global K;

%b
f = @(x) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - 11 * x .* exp(-x);
f_prim = @(x) (21*(x.^2 + x + 0.05).^7) ./ (3*x + 1).^8 - (7*(2*x + 1).*(x.^2 + x + 0.05).^6) ./ (3*x + 1).^7 - 11*exp(-x) + 11*x.*exp(-x) + 73;

%c
function root = find_root_newton(x_0)
    %Import mathematical functions
    global f;
    global f_prim;
    global K;

    x_root = x_0;
    error_prev = 1;

    for i = 1:100
        x_delta = f(x_root) / f_prim(x_root);
        
        %Error and convergence testing
        error = abs(x_delta);
        x_root = x_root - x_delta;
        relative_error = error / x_root;

        K = error / (error_prev^2);
        error_prev = error;
        disp(["X:", x_root, "   K:" , K,]);

        %Exit condition
        if relative_error <= 10^(-8)
            disp(["iterations: ", i]);
            break;
        end
    end
    root = [x_root, relative_error];
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
roots = [find_root_newton(0); find_root_newton(7)];
plot(roots, f(roots), "o");

%Plot first root
subplot(2, 2, 2)
x = linspace(roots(1)-0.01, roots(1)+0.01, 500);
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
%Kvadratisk konvergens definition: lim då (n -> inf); abs(error) / (abs(prev_error)^2) = K


%E
%Bestäm konvergenskonstanten.

%x_root = roots(2);

%x_delta = f(x_root) / f_prim(x_root);
%error1 = abs(x_delta);
%error2 = abs(f(x_root - x_delta) / f_prim(x_root - x_delta));
%K = error2 / (error1^2);

K


%F och G. Visa rötter och att konvergens var OK. Konvergens ser man i att K-värdet hålls ungefär samma mellan iterationerna.
disp("Root:                 Relative error:");
disp(roots);