function ret = own_ode45(odefun, x_span, y0, number_of_steps, tolerance)
    n = number_of_steps;
    num_equations = length(y0);
    prev_last_y = y0;

    for iteration = 1:20
        h = (x_span(2) - x_span(1)) / n;
        xx = x_span(1):h:x_span(2);
        yy = zeros(num_equations, length(xx));
        yy(1) = y0;

        %Do Runge-Kutta solve
        for i = 1:(n-1)
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

        n = n * 2;
        prev_last_y(yy(:,end))
    end
end

function ret = own_polyfit(x_coords, y_coords, grade)

end

function ret = own_polyval(polynom, x_query)

end

function ret = own_spline(x_coords, y_coords, x_query)

end