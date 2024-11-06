clear all; clc; close all;
import Rocket.*;
f1=figure;


rocket1=Rocket(0,0,0,0,45,20,3,0.01,10);

%Solve and plot big rocket's trajectory
rocket1=rocket1.solve_trajectory(20);
plot(rocket1.x_pos,rocket1.y_pos,'-o')
hold on 
plot([0,70],[0,0])

%Find and plot the max height point of large rocket
[max_height_x, max_height_y, max_height_index] = rocket1.get_highest_point(); %SVAR på uppgift a
hold on;
plot(max_height_x, max_height_y, "o");

%Find and plot landing point of large rocket
land_point = rocket1.get_land_point;
hold on
plot(land_point(1),land_point(2),'o')
disp(land_point)

%Solve small rocket best fire point/time.
n=10;
little_rockets= createArray(1,n,"Rocket");

%Interpolate start conditions for small rocket.
tt=linspace(0,rocket1.t_values(max_height_index+1),n);
little_rocket_y_start=spline(rocket1.t_values,rocket1.y_pos,tt);
little_rocket_x_start=spline(rocket1.t_values,rocket1.x_pos,tt);
little_rocket_y_vel=spline(rocket1.t_values,rocket1.y_vel,tt);
little_rocket_x_vel=spline(rocket1.t_values,rocket1.x_vel,tt);

%Solve and plot small rocket trajectories
f2=figure;
plot(rocket1.x_pos, rocket1.y_pos, ":")
for i = 1:length(tt)
    start_x_pos=little_rocket_x_start(i);
    start_y_pos=little_rocket_y_start(i);
    start_x_vel=little_rocket_x_vel(i);
    start_y_vel=little_rocket_y_vel(i);
    %Rocket(start_x_pos,start_y_pos,start_x_vel,start_y_vel,12,5,4,0.001,3);
    little_rockets(i)=Rocket(start_x_pos,start_y_pos,start_x_vel,start_y_vel,12,5,4,0.001,3).solve_trajectory(40);
    hold on
    plot(little_rockets(i).x_pos,little_rockets(i).y_pos,'o-')
    %hold on
end

hold on
plot([0,300], [0,0]);

%Get landing points of small rocket
land_points=createArray(2,n);
for i = 1:n
    land_points(:,i)=little_rockets(i).get_land_point();
end

hold on
plot(land_points(1,:), land_points(2,:), "o");