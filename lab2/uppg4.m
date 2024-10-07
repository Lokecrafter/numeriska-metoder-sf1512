close all; clear all; clc;


points = [19, 16, 5; 
           5, 19, 16];

circle_origin = [0; 0];
circle_radius = 1;


function ret = plot_cirlce(origin, radius);
    tt = linspace(0, 2 * pi, 100);

    xx = origin(1) + cos(tt) * radius;
    yy = origin(2) + sin(tt) * radius;

    plot(xx, yy, "-r");
    hold on;
    plot(origin(1), origin(2), "or");
end

plot(points(1,:), points(2,:), "o");
hold on;

circle_origin(1) = mean(points(1,:));
circle_origin(2) = mean(points(2,:));

distances_x_y_squared = (points - circle_origin).^2;
distances = sqrt(distances_x_y_squared(1,:) + distances_x_y_squared(2,:));
circle_radius = mean(distances);


prev_circle_origin = circle_origin;
for i = 1:100
    new_circle_origin = circle_origin;
    new_radius = 0;

    for point_index = 1:length(points)
        point_direction = points(:,point_index) - circle_origin;
        distance_to_point = norm(point_direction);
        point_direction = point_direction / distance_to_point;


        new_circle_origin = new_circle_origin + point_direction * (distance_to_point - circle_radius) * 0.5;
        new_radius = new_radius + distance_to_point;
    end

    circle_radius = new_radius / 3;
    circle_origin = new_circle_origin;

    E_trunk = norm(circle_origin - prev_circle_origin);
    prev_circle_origin = circle_origin;

    disp("E_trunk: " + E_trunk);
end


plot_cirlce(circle_origin ,circle_radius)
pbaspect([1 1 1]) %Sets square aspect ratio