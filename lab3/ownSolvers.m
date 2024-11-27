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
    h = zeros(1,length(x_coords)-1);
    yh = zeros(1,length(y_coords)-1);
    c = zeros(1,length(y_coords)-1);
    for i = 1:length(h)
        h(i) = x_coords(i+1)-x_coords(i);
        yh(i) = y_coords(i+1)-y_coords(i);
        c(i) = yh(1)/h(i);
        
    end
    A = zeros(length(x_coords),length(x_coords));
    b = zeros(length(x_coords),1);
    for row = 2:length(x_coords)-1
        A(row,row)=2*(h(row)+h(row-1));
        A(row,row+1)=h(row-1);
        A(row,row-1)=h(row);
        b(row) = yh(i) * h(i-1)/h(i) + yh(i-1) * h(i)/h(i-1);
    end
    A(1,1)=2*h(1);
    A(1,2)=h(1);
    A(end,end)=2*h(end);
    A(end,end-1)=h(end);
    b(1) = yh(1);
    b(end) = yh(end);

    disp(h)
    disp(A)
    disp(b)

    k = A\(3*b);
    y_query = zeros(size(x_query));
    for j = 1:length(x_query)
        x = x_query(j);
        for i = 1:length(x_coords)-1
            is_inside_interval = discretize(x, [x_coords(i),x_coords(i+1)])==1;
            if is_inside_interval 
                disp(i)
                break
            end
        end
        y_query(j) = y_coords(i) + c(i)*(x-x_coords(i)) + (x-x_coords(i))*(x-x_coords(i+1))*((k(i+1)-c(i))*(x-x_coords(i))+(k(i)-c(i))*(x-x_coords(i+1)))/(h(i)*h(i));
    end
    ret = y_query;
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
x_data = [-1, 0, 1, 2];
y_data = [1, 2, 1, 3];
plot(x_data, y_data, "o-")

% p = polyfit(x_data, y_data, 2);
% disp(p)
% hold on
% plot(xx, polyval(p,xx), "o");

%p = own_polyfit(x_data, y_data, 2);
%disp(p)

%hold on
%plot(xx, polyval(p,xx), "o");
%hold on
%plot(xx, own_polyval(p,xx), ":");
hold on
yy = own_spline(x_data,y_data,xx);
plot(xx,yy,'o:')
hold on
yy = spline(x_data,y_data,xx);
plot(xx,yy,'o-')

