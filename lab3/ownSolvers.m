disp("Hello!");

function ret = own_ode45(odefun, x_span, y0, number_of_steps, tolerance)
    n = number_of_steps;
    num_equations = length(y0);
    prev_last_y = y0;

    for iteration = 1:20
        h = (x_span(2) - x_span(1)) / n;
        xx = x_span(1):h:x_span(2);
        yy = zeros(num_equations, length(xx));
        yy(:,1) = y0;

        %Do Runge-Kutta solve
        for i = 1:n
            k1 = odefun(xx(i), yy(i));
            k2 = odefun(xx(i) + h * 0.5, yy(i) + h * k1 * 0.5);
            k3 = odefun(xx(i) + h * 0.5, yy(i) + h * k2 * 0.5);
            k4 = odefun(xx(i) + h,       yy(i) + h * k3      );

            K = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            yy(:,i+1) = yy(:,i) + K*h;
        end

        y_diff = yy(:,end) - prev_last_y;
        E_trunk = norm(y_diff);

        %Exit condition
        if iteration >= 2
            if E_trunk < tolerance
                break
            end
        end

        prev_last_y = yy(:,end);
        n = n * 2;
    end
    ret.x = xx;
    ret.y = yy;
    ret.E_trunk = E_trunk;
    ret.iterations = iteration;
end

function ret = own_polyfit(x_coords, y_coords, grade)
    A = ones(length(x_coords), grade+1);
    for i = 1:length(x_coords)
        A(i,1:end-1) = x_coords(i);
        for j = 1:grade
            A(i,j) = A(i,j) .^ (grade + 1 - j);
        end
    end

    b = y_coords';
    c = A\b;

    ret = c';
end

function ret = own_polyval(polynom, x_query)
    y_result = zeros(size(x_query));
    grade = length(polynom) - 1;

    for i = 1:grade+1
        y_result = y_result + polynom(i) .* (x_query .^ (grade + 1 - i));
    end

    ret = y_result;
end

function ret = own_spline(x_coords, y_coords, x_query)

end



% odefun = @(t, y) 0 * y + t;
% x_span = [0, 1];
% y0 = [0; 0];

% result = own_ode45(odefun, x_span, y0, 10, 1e-2);
% plot(result.x, result.y);
% hold on
% disp("E_trunk: " + result.E_trunk + "   Iterations: " + result.iterations);

% result = ode45(odefun, x_span, y0);
% plot(result.x, result.y);


xx = linspace(-1, 2, 10);
x_data = [-1, 2, 1, 2];
y_data = [1, 2, 1, 3];
plot(x_data, y_data, "o-")

% p = polyfit(x_data, y_data, 2);
% disp(p)
% hold on
% plot(xx, polyval(p,xx), "o");

p = own_polyfit(x_data, y_data, 2);
disp(p)

hold on
plot(xx, polyval(p,xx), "o");
hold on
plot(xx, own_polyval(p,xx), ":");