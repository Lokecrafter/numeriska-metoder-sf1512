close all, clear all, clc;

%Minsta positiva roten är nästan noll.

%Tar bort tal flera gånger så får man roten:
%7.37

f = @(x) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - 11 * x .* exp(-x);
f_prim = @(x) (21*(x.^2 + x + 0.05).^7) ./ (3*x + 1).^8 - (7*(2*x + 1).*(x.^2 + x + 0.05).^6) ./ (3*x + 1).^7 - 11*exp(-x) + 11*x.*exp(-x) + 73;


function root = find_root_newton(x_0)
    f = @(x) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - 11 * x .* exp(-x);
    f_prim = @(x) (21*(x.^2 + x + 0.05).^7) ./ (3*x + 1).^8 - (7*(2*x + 1).*(x.^2 + x + 0.05).^6) ./ (3*x + 1).^7 - 11*exp(-x) + 11*x.*exp(-x) + 73;

    x_prev = x_0;
    for i = 1:100
        x_delta = f(x_prev) / f_prim(x_prev);
        x_prev = x_prev - x_delta;
    end
    root = x_prev;
end




x = linspace(-1, 8, 500);

%Plot function
plot(x, f(x));
ylim([-500, 500]);
hold on;

%Plot derivative
plot(x, f_prim(x));
hold on;

%Plot zero reference line
plot(x, zeros(length(x)));
hold on;

%Calc and plot roots
roots = [find_root_newton(0), find_root_newton(7)]
plot(roots, f(roots), "o");