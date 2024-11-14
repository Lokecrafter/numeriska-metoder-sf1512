clear all; clc; close all;
import Rocket.*;
f1=figure;
pbaspect([1 1 1])

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
number_of_small_rockets=20;
little_rockets= createArray(1,number_of_small_rockets,"Rocket");

t_interval=[0,rocket1.t_values(max_height_index+1)];

function ret=find_three_best_candidates(t_values,x_values)
    [~,index1]=max(x_values);
    % x_values(index1)=0;
    % [~,index2]=max(x_values);
    % x_values(index2)=0;
    % [~,index3]=max(x_values);
    % x_values(index3)=0;
    %ret=[t_values(index1),t_values(index2),t_values(index3)];
    ret=[t_values(index1),t_values(index1-1),t_values(index1 +1)];
end

f2=figure;
for iteration =1:4

    %Interpolate start conditions for small rocket.
    small_rocket_start_times=linspace(t_interval(1),t_interval(2),number_of_small_rockets);
    little_rocket_y_start=spline(rocket1.t_values,rocket1.y_pos,small_rocket_start_times);
    little_rocket_x_start=spline(rocket1.t_values,rocket1.x_pos,small_rocket_start_times);
    little_rocket_y_vel=spline(rocket1.t_values,rocket1.y_vel,small_rocket_start_times);
    little_rocket_x_vel=spline(rocket1.t_values,rocket1.x_vel,small_rocket_start_times);

    %Solve and plot small rocket trajectories
    %f2=figure;
    %plot(rocket1.x_pos, rocket1.y_pos, ":")
    for i = 1:number_of_small_rockets
        start_x_pos=little_rocket_x_start(i);
        start_y_pos=little_rocket_y_start(i);
        start_x_vel=little_rocket_x_vel(i);
        start_y_vel=little_rocket_y_vel(i);
        %Rocket(start_x_pos,start_y_pos,start_x_vel,start_y_vel,12,5,4,0.001,3);
        little_rockets(i)=Rocket(start_x_pos,start_y_pos,start_x_vel,start_y_vel,12,5,4,0.001,3).solve_trajectory(40);
        hold on
       %plot(little_rockets(i).x_pos,little_rockets(i).y_pos,'o-')
        %hold on
    end

    
    %Get landing points of small rocket
    land_points=createArray(2,number_of_small_rockets);
    for i = 1:number_of_small_rockets
        land_points(:,i)=little_rockets(i).get_land_point();
    end

    hold on
    plot(small_rocket_start_times,land_points(1,:),'-o')

    %disp(max(land_points(1,:)))

    [x_value,index]=max(land_points(1,:));
    disp(x_value)
    %t_interval(1)=small_rocket_start_times(index-1);
    %t_interval(2)=small_rocket_start_times(index+1);
    best_times=find_three_best_candidates(small_rocket_start_times,land_points(1,:));
    t_interval(1)=min(best_times);
    t_interval(2)=max(best_times);

    
end

[~, max_x_distance_index] = max(land_points(1,:));
max_x_distance_small_rocket_polynom = polyfit(small_rocket_start_times(max_x_distance_index-1:max_x_distance_index+1),land_points(1,max_x_distance_index-1:max_x_distance_index+1),2);

tt=linspace(2.94, 3.06, 1000);
hold on
plot(tt, polyval(max_x_distance_small_rocket_polynom, tt), "-.");

p = max_x_distance_small_rocket_polynom;
time_to_fire_small_rocket = -p(2)/(2*p(1));
disp("Time to fire small rocket: " + time_to_fire_small_rocket);

hold on
plot([time_to_fire_small_rocket, time_to_fire_small_rocket], [-200, 300], "-")

hold on
plot(rocket1.t_values,rocket1.y_pos,'-o')

f3=figure;
pbaspect([1 1 1])
hold on
plot([0,300], [0,0]);

hold on
plot(land_points(1,index), land_points(2,index), "o")

hold on
plot(little_rockets(index).x_pos,little_rockets(index).y_pos,'o-')

hold on 
plot(rocket1.x_pos, rocket1.y_pos, ":")
