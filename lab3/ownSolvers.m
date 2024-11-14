function ret = own_ode45(odefun, x_span, y0, number_of_steps, tolerance)
    n = 10;
    h = (x_span(2) - x_span(1)) / n;
    xx = x_span(1):h:x_span(2);
    yy = zeros(length(y0), length(xx));
    yy(1) = y0;

    for i = 1:(n-1)
        k1 = odefun(xx(i), yy(i));
        k2 = odefun(xx(i) + h * 0.5, yy(i) + h * k1 * 0.5);
        k3 = odefun(xx(i) + h * 0.5, yy(i) + h * k2 * 0.5);
        k4 = odefun(xx(i) + h,       yy(i) + h * k3      );

        K = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        y(i+1) = y(i) + K*h;
    end
end

function ret = own_polyfit(x_coords, y_coords, grade)

end

function ret = own_polyval(polynom, x_query)

end

function ret = own_spline(x_coords, y_coords, x_query)

end