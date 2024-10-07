%upg 4 lab 2

close all; clear all; clc;

global points

points = [19, 16, 5; 
           5, 19, 16]; % x and y given as columns

circle_origin = [0; 0];
circle_radius = 1;

tolerance = 1e-5;

% function which returns the jacobian matrix
function ret = get_jacobian(X, Y, R);
    global points
    
    number_of_equations = length(points); % length of matrix gives the length of the longest dimension (row vs column).
    number_of_variables = 3;
    Jacobian_matrix = zeros([number_of_equations, number_of_variables]);

    for equation_index = 1:number_of_equations
        Jacobian_matrix(equation_index, 1) = - 2 * (points(1, equation_index) - X); %df/dX part. der.
        Jacobian_matrix(equation_index, 2) = - 2 * (points(2, equation_index) - Y); %df/dY part. der.
        Jacobian_matrix(equation_index, 3) = - 2 * R; %df/dR part. der.
    end

    ret = Jacobian_matrix;
end

% plot circle funct.
function ret = plot_cirlce(origin, radius);
    tt = linspace(0, 2 * pi, 100);

    xx = origin(1) + cos(tt) * radius;
    yy = origin(2) + sin(tt) * radius;

    plot(xx, yy, "-r");
    hold on;
    plot(origin(1), origin(2), "or");
end

% also plots circle
plot(points(1,:), points(2,:), "o");
hold on;

start_guess_NR = [0,0,1]';

% Right side of equation f(X,Y,R) = 
function ret = get_function_values(X, Y, R);
    global points
    
    number_of_equations = length(points); % length of matrix gives the length of the longest dimension (row vs column).
    function_values = zeros([number_of_equations, 1]);

    for i = 1:number_of_equations
        function_values(i) = (points(1, i) - X)^2 + (points(2, i) - Y)^2 - R^2;
    end

    ret = function_values;

end

%-----------------------------------------------
%
% NR-method for equation-system
%
%-----------------------------------------------

current_guess = start_guess_NR;
prev_guess = 0;

for i = 1:10

    curr_X = current_guess(1);
    curr_Y = current_guess(2);
    curr_R = current_guess(3);

    jacobian_matrix = get_jacobian(curr_X, curr_Y, curr_R);
    function_values = - get_function_values(curr_X, curr_Y, curr_R);

    t = jacobian_matrix\function_values;

    current_guess = current_guess - t;

    E_trunk = abs(current_guess - prev_guess);

    % exit cond.
    if max(E_trunk) < tolerance % max value of error vector
        if i >=2
            break
        end
    end

    prev_guess = current_guess;

end

circle_origin(1) = current_guess(1);
circle_origin(2) = current_guess(2);
circle_radius = current_guess(3);


% square aspect for circle plot
plot_cirlce(circle_origin, circle_radius)
pbaspect([1 1 1]) %Sets square aspect ratio