classdef Solvers
    properties
        use_built_in
    end
    methods
        function obj=Solvers(use_matlab_functions)
            if nargin == 1
                obj.use_built_in = use_matlab_functions;
            end
        end

        function ret = solve_ode45(obj, odefun, x_span, y0, number_of_steps, tolerance)
            if obj.use_built_in
                options1 = odeset('RelTol',tolerance);
                result1 = ode45(odefun, x_span, y0, options1);
                options2 = odeset('RelTol',tolerance * 0.5);
                result2 = ode45(odefun, x_span, y0, options2);

                E_trunk = abs(result2.y(1, length(result2.y) - 1) - result1.y(1, length(result1.y) - 1));

                ret = result2;
                ret.E_trunk = E_trunk;
                return;
            end
            n = number_of_steps;
            num_equations = length(y0);
            prev_last_y = y0;
            last_E_trunk = 0;

            for iteration = 1:7
                h = (x_span(2) - x_span(1)) / n;
                xx = x_span(1):h:x_span(2);
                yy = zeros(num_equations, length(xx));
                yy(:,1) = y0;

                %Do Runge-Kutta solve
                for i = 1:n
                    k1 = odefun(xx(i), yy(:,i));
                    k2 = odefun(xx(i) + h * 0.5, yy(:,i) + h * k1 * 0.5);
                    k3 = odefun(xx(i) + h * 0.5, yy(:,i) + h * k2 * 0.5);
                    k4 = odefun(xx(i) + h,       yy(:,i) + h * k3      );

                    K = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                    yy(:,i+1) = yy(:,i) + K*h;
                end

                y_diff = yy(:,end) - prev_last_y;
                E_trunk = norm(y_diff);

                disp("E_trunk: " + E_trunk + "   Last E_trunk: " + last_E_trunk + "    K: " + (last_E_trunk / (E_trunk.^4)));

                %Exit condition
                if iteration >= 2
                    if E_trunk < tolerance
                        break
                    end
                end

                prev_last_y = yy(:,end);
                last_E_trunk = E_trunk;
                n = n * 2;
            end
            ret.x = xx;
            ret.y = yy;
            ret.E_trunk = E_trunk;
            ret.iterations = iteration;
        end

        function ret = solve_polyfit(obj, x_coords, y_coords, grade)
            if obj.use_built_in
                ret = polyfit(x_coords, y_coords, grade);
                return;
            end

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

        function ret = solve_polyval(obj, polynom, x_query)
            if obj.use_built_in
                ret = polyval(polynom, x_query);
                return;
            end

            y_result = zeros(size(x_query));
            grade = length(polynom) - 1;

            for i = 1:grade+1
                y_result = y_result + polynom(i) .* (x_query .^ (grade + 1 - i));
            end

            ret = y_result;
        end

        function ret = solve_interpolate(obj, x_coords, y_coords, x_query)
            if obj.use_built_in
                ret = interp1(x_coords, y_coords, x_query);
                return;
            end

            yy = zeros(size(x_query));

            for query_i = 1:length(x_query)
                for i = 1:length(x_coords)
                    if x_query(query_i) < x_coords(i) %Will automatically make i=length if none is found
                        break
                    end
                end

                x1 = x_coords(i - 1);
                x2 = x_coords(i);
                y1 = y_coords(i - 1);
                y2 = y_coords(i);

                %Linear interpolation
                yy(query_i) = ((x_query(query_i) - x1) / (x2 - x1)) * (y2 - y1) + y1;
            end

            ret = yy;
        end

        function ret = solve_spline(obj, x_coords, y_coords, x_query)
            xh = zeros(1, length(x_coords)-1);
            yh = zeros(1, length(y_coords)-1);
            c =  zeros(1, length(y_coords)-1);

            %Calc constants 
            for i = 1:length(xh)
                xh(i) = x_coords(i+1) - x_coords(i);
                yh(i) = y_coords(i+1) - y_coords(i);
                c(i) = yh(i) / xh(i);
            end

            %Construct splines matrix and vector
            A = zeros(length(x_coords));
            b = zeros(length(x_coords),1);
            for row = 2:length(x_coords)-1
                A(row,row-1) = xh(row);
                A(row,row) = 2 * (xh(row) + xh(row-1));
                A(row,row+1) = xh(row-1);
                b(row) = yh(i) * xh(i-1)/xh(i) + yh(i-1) * xh(i)/xh(i-1);
            end

            A(1,1) = 2*xh(1);
            A(1,2) = xh(1);
            A(end,end) = 2*xh(end);
            A(end,end-1) = xh(end);
            b(1) = yh(1);
            b(end) = yh(end);

            % disp(xh)
            % disp(A)
            % disp(b)

            %K-values in Hermite's
            k = A\(3*b);

            y_query = zeros(size(x_query));
            for j = 1:length(x_query)
                x = x_query(j);
                for i = 1:length(x_coords)-1
                    is_inside_interval = discretize(x, [x_coords(i),x_coords(i+1)])==1;
                    if is_inside_interval 
                        break
                    end
                end
                y_query(j) = y_coords(i) + c(i)*(x-x_coords(i)) + (x-x_coords(i)) * (x-x_coords(i+1)) * ((k(i+1) - c(i)) * (x - x_coords(i)) + (k(i) - c(i)) * (x - x_coords(i+1))) / (xh(i)*xh(i));
            end
            ret = y_query;
        end
    end
end