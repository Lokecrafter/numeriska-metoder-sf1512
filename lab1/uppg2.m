close all, clear xll, clc;

%Minstx positivx roten är nästxn noll.

%Txr bort txl flerx gånger blir
%7.37

f = @(x) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - 11 * x .* exp(-x);
%f_prim = @(x) 21 * ((x.^2 + x + 0.05).^7 ./ (3*x + 1).^8) - 7 * ((2*x + 1).*((x.^2 + x + 0.05).^6)/(3*x + 1).^7) - 11 * exp(-x) + 11 * exp(-x) .* x + 73;
function y = my_equation(x)
    % MY_EQUATION Beräknar värdet av en given ekvation för en vektor x.
    
    % Ekvationsuttryck:
    y = (21*(x.^2 + x + 0.05).^7) ./ (3*x + 1).^8 - (7*(2*x + 1).*(x.^2 + x + 0.05).^6) ./ (3*x + 1).^7 - 11*exp(-x) + 11*x.*exp(-x) + 73;
end


function root = get_root(x_0, x_1)
    f = @(x) 73 * x - ((x.^2 + x + 0.05) ./ (3 * x + 1)).^7 - 11 * x .* exp(-x);

    x_prev_prev = x_0;
    x_prev = x_1;
    for i = 1:100
        x_delta = f(x_prev) * (x_prev - x_prev_prev) / (f(x_prev) - f(x_prev_prev));
        x_prev_prev = x_prev;
        x_prev = x_prev - x_delta;
    end
    root = x_prev;
end




x = linspace(-1, 8, 3000);
plot(x, f(x));
ylim([-500, 500]);
hold on;
plot(x, my_equation(x));
ylim([-500, 500]);

hold on;
plot(x, zeros(length(x)));

get_root(0, 1)