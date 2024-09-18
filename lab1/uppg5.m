close all, clear all, clc;



global f;
global f_prim;
global K;
global modify_factor

modify_factor = 1;

%b
f = @(x, q) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - q * 11 * x .* exp(-x);
f_prim = @(x, q) (21*(x.^2 + x + 0.05).^7) ./ (3*x + 1).^8 - (7*(2*x + 1).*(x.^2 + x + 0.05).^6) ./ (3*x + 1).^7 - q*11*exp(-x) + q*11*x.*exp(-x) + 73;

%c
function root = find_root_newton(x_0)
    %Import mathematical functions
    global f;
    global f_prim;
    global K;
    global modify_factor

    modify_factor
    x_root = x_0;
    error_prev = 1;
    error = 1;
    
    for i = 1:100
        x_delta = f(x_root, modify_factor) / f_prim(x_root, modify_factor);
        
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
    root = [x_root, error];
end



%B plotting function--------------------------------------------------------------------------------------------
x = linspace(-1, 8, 500);

%Plot function
plot(x, f(x, modify_factor));
title("f(x)");
grid on;
ylim([-500, 500]);
hold on;
%Plot zero reference line
plot(x, zeros(length(x)));
hold on;

%C Find roots and plot them-----------------------------------------------------------------------------------------
modify_factor = 1;
roots_original = [find_root_newton(0); find_root_newton(7)];
plot(roots_original, f(roots_original, modify_factor), "o");
hold on;

modify_factor = 1.5;
roots_risen = [find_root_newton(0); find_root_newton(7)];
plot(roots_risen, f(roots_risen, modify_factor), "o");
hold on;

modify_factor = 0.75;
roots_lowered = [find_root_newton(0); find_root_newton(7)];
plot(roots_lowered, f(roots_lowered, modify_factor), "o");




%F och G. Visa rötter och att konvergens var OK. Konvergens ser man i att K-värdet hålls ungefär samma mellan iterationerna.
disp("Root:                 Absolute error:");
disp(roots_original);
disp(roots_risen);
disp(roots_lowered);

disp("Below is in %");

disp("Constant raised by 3%:");
disp([((roots_risen ./ roots_original) - 1) * 100]);
disp("Constant lowered by 3%:");
disp([((roots_lowered ./ roots_original) - 1) * 100]);