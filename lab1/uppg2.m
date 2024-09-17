close all, clear all, clc;


%a
%Minsta positiva roten är nästan noll.

%Största rot
%Tar bort tal flera gånger så får man roten:
%7.37
global f;
global f_prim;

%b
f = @(x) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - 11 * x .* exp(-x);
f_prim = @(x) (21*(x.^2 + x + 0.05).^7) ./ (3*x + 1).^8 - (7*(2*x + 1).*(x.^2 + x + 0.05).^6) ./ (3*x + 1).^7 - 11*exp(-x) + 11*x.*exp(-x) + 73;

%c
function root = find_root_newton(x_0)
    global f;
    global f_prim;

    x_prev = x_0;
    error_prev = 1;
    for i = 1:100
        x_delta = f(x_prev) / f_prim(x_prev);
        error = x_delta;

        K = error / (error_prev^2)
        if error <= 10^(-8)
            i
            break;
        end

        x_prev = x_prev - x_delta;

        error_prev = error;
    end
    root = x_prev;
end




x = linspace(-1, 8, 500);

%Plot function
subplot(2, 2, 1)
plot(x, f(x));
ylim([-500, 500]);
hold on;

%Plot zero reference line
plot(x, zeros(length(x)));
hold on;

%Find roots
roots = [find_root_newton(0), find_root_newton(7)];
plot(roots, f(roots), "o");

%Plot first root
subplot(2, 2, 2)
x = linspace(roots(1)-0.1, roots(1)+0.1, 500);
plot(x, f(x));
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
hold on;
plot(roots(2), f(roots(2)), "o");
hold on;
%Plot zero reference line
plot(x, zeros(length(x)));
hold on;


%D
%Kvadratisk konvergens definition: lim då (n -> inf); abs(error) / (abs(prev_error)^2) = K


%E
%Bestäm konvergenskonstanten.
K = 0;
x_prev = roots(2
);
error_prev = 1;
for i = 1:100
    x_delta = f(x_prev) / f_prim(x_prev);
    error = x_delta;

    K = error / (error_prev^2)
    if error <= 10^(-8)
        i
        break;
    end

    x_prev = x_prev - x_delta;

    error_prev = error;
end
root = x_prev;

K

%1.3213e-12