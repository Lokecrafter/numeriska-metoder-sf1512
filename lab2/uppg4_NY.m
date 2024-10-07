%upg 4 lab 2

close all; clear all; clc;

global points

points = [19, 16, 5; 
          5, 19, 16]; % x and y given as columns
% points = [19, 16, 5, 24, 2; 
%            5, 19, 16, 2, 20]; % x and y given as columns

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
function ret = plot_cirlce(origin, radius, plot_settings);
    tt = linspace(0, 2 * pi, 100);

    xx = origin(1) + cos(tt) * radius;
    yy = origin(2) + sin(tt) * radius;

    plot(xx, yy, "-" + plot_settings);
    hold on;
    plot(origin(1), origin(2), "or");
end

% also plots circle
plot(points(1,:), points(2,:), "o");
hold on;

start_guess_NR = [0,0,8]';

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

for i = 1:20
    curr_X = current_guess(1);
    curr_Y = current_guess(2);
    curr_R = abs(current_guess(3));

    jacobian_matrix = get_jacobian(curr_X, curr_Y, curr_R);
    function_values = get_function_values(curr_X, curr_Y, curr_R);

    %Multiply first equation by five
    %jacobian_matrix(1,:) = jacobian_matrix(1,:) .* 5;
    %function_values(1) = function_values(1) .*5;

    t = jacobian_matrix\function_values;

    current_guess = current_guess - t;

    E_trunk = abs(current_guess - prev_guess);



    disp("X: " + current_guess(1) + "    Y: " + current_guess(2) + "    R: " + current_guess(3) + "    E_trunk: " + max(E_trunk));
    % exit cond.
    if max(E_trunk) < tolerance % max value of error vector
        if i >=2
            disp("Iterations: " + i);
            break
        end
    end

    prev_guess = current_guess;
end

circle_origin(1) = current_guess(1);
circle_origin(2) = current_guess(2);
circle_radius = current_guess(3);


% square aspect for circle plot
plot_cirlce(circle_origin, circle_radius, "r")
pbaspect([1 1 1]) %Sets square aspect ratio






A = zeros([length(points), 3]);
A(:,1) = ones([1, length(points)]);
A(:,2) = points(1,:)';
A(:,3) = points(2,:)';

b = (points(1,:).^2 + points(2,:).^2)';

c = A\b;

X = 0.5 * c(2);
Y = 0.5 * c(3);
R = sqrt(c(1) + 0.25 * c(2) * c(2) + 0.25 * c(3) * c(3));


hold on;
plot_cirlce([X,Y], R, "g");


