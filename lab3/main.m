clear all; clc; close all; format long;
import Rocket.*;
import Solvers.*;

%Inputs
use_matlab_functions_for_solvers = true;
tolerance = 0.01;

rocket1_mass = 20;
rocket1_fuel = 45;
rocket1_burn_time = 3;
rocket1_burn_force = 10;
rocket1_air_resistance = 0.01;
rocket1_trajectory_end_time = 30;

rocket2_mass = 5;
rocket2_fuel = 12;
rocket2_burn_time = 4;
rocket2_burn_force = 3;
rocket2_air_resistance = 0.001;
rocket2_trajectory_end_time = 50;
number_of_small_rockets=20; %Used for calculating best launch point

%Numerical methods
own_solver = Solvers(use_matlab_functions_for_solvers);




%Main program
f1=figure;
disp("Program start" + newline);

disp("Calculating max height and landing point for large rocket...");
%Initialize big rocket
angle = deg2rad(80);
start_vel = 20 * [cos(angle), sin(angle)];
rocket1=Rocket(0,0,start_vel(1),start_vel(2),rocket1_fuel,rocket1_mass,rocket1_burn_time,rocket1_air_resistance,rocket1_burn_force);

%Solve and plot big rocket's trajectory
rocket1=rocket1.solve_trajectory(rocket1_trajectory_end_time,tolerance);
hold on;    plot(rocket1.x_pos,rocket1.y_pos,'-*')
hold on;    plot([0,32],[0,0])
disp("Large rocket global trajectory E_trunk: " + sprintf('%0.10f', rocket1.get_trajectory_E_trunk()));

%Find and plot the max height point of large rocket
rocket1_max_height = rocket1.get_highest_point(); %SVAR på uppgift a
hold on;    plot(rocket1_max_height.point(1), rocket1_max_height.point(2), "o",'MarkerSize',10);
disp("Large rocket max height.      x: " + sprintf('%0.6f', rocket1_max_height.point(1)) + "   y: " + sprintf('%0.6f', rocket1_max_height.point(2)) + "   E_trunk: " + sprintf('%0.10f', rocket1_max_height.E_trunk) + "   Global E_trunk: " + sprintf('%0.10f', rocket1_max_height.glob_E_trunk))

%Find and plot landing point of large rocket
land_point = rocket1.get_land_point();
hold on;    plot(land_point.point(1),land_point.point(2),'o','MarkerSize',10)
disp("Large rocket landing point.   x: " + sprintf('%0.6f', land_point.point(1)) + "   y: " + sprintf('%0.6f', land_point.point(2))+ "   E_trunk: " + sprintf('%0.10f', land_point.E_trunk) + "   Global E_trunk: " + sprintf('%0.10f', land_point.glob_E_trunk))

%Make plots pretty
title("Large rocket trajectory");
xlabel("x position");
ylabel("y position");
legend(["Trajectory", "Ground", "Max height.   x: " + rocket1_max_height.point(1) + ",   y:" + rocket1_max_height.point(2), "Land point.   x: " + land_point.point(1) + ",   y:" + land_point.point(2)], 'Location', 'northwest');


function ret=find_three_best_candidates(t_values,x_values)
    [~,index1]=max(x_values);
    if index1 + 1 > length(t_values)
        ret.t_values=[t_values(index1 - 2),t_values(index1-1),t_values(index1)];
        ret.indexes = [index1-2, index1-1, index1];
        return;
    else if index1 - 1 < 1
        ret.t_values=[t_values(index1),t_values(index1+1),t_values(index1+2)];
        ret.indexes = [index1, index1+1, index1+2];
        return;
        end
    end
    
    ret.t_values = [t_values(index1),t_values(index1-1),t_values(index1+1)];
    ret.indexes = [index1, index1-1, index1+1];
end
%{
    %Plot big rocket solutions against time
    % figureROcket1Interp = figure;
    %     tspan1=[0,rocket1.burn_time];
    %     tspan2=[rocket1.burn_time,rocket1.burn_time+2];
    %     tspan3=[rocket1.burn_time+2,30];
    
    %     val_amount = 1000;
    %     tt = [linspace(tspan1(1), tspan1(2), val_amount), linspace(tspan2(1), tspan2(2), val_amount), linspace(tspan3(1), tspan3(2), val_amount)];
    
    %     subplot(2, 2, 1);
    %     hold on;    plot(rocket1.t_values, rocket1.y_pos, 'o');
    %     y_start=own_solver.solve_interpolate(rocket1.t_values,rocket1.y_pos,tt);
    %     % y_start=own_solver.solve_spline(rocket1.t_values,rocket1.y_pos,tt);
    %     hold on;    plot(tt, y_start, ':');
    %     y_start=spline(rocket1.t_values,rocket1.y_pos,tt);
    %     plot(tt, y_start, '-');

    %     subplot(2, 2, 2);
    %     hold on;    plot(rocket1.t_values, rocket1.x_pos, 'o');
    %     x_start=own_solver.solve_interpolate(rocket1.t_values,rocket1.x_pos,tt);
    %     % x_start=own_solver.solve_spline(rocket1.t_values,rocket1.x_pos,tt);
    %     hold on;    plot(tt, x_start, ':');
    %     x_start=spline(rocket1.t_values,rocket1.x_pos,tt);
    %     plot(tt, x_start, '-');

    %     subplot(2, 2, 3);
    %     hold on;    plot(rocket1.t_values, rocket1.y_vel, 'o');
    %     y_vel=own_solver.solve_interpolate(rocket1.t_values,rocket1.y_vel,tt);
    %     % y_vel=own_solver.solve_spline(rocket1.t_values,rocket1.y_vel,tt);
    %     hold on;    plot(tt, y_vel, ':');
    %     y_vel=spline(rocket1.t_values,rocket1.y_vel,tt);
    %     hold on;    plot(tt, y_vel, '-');

    %     subplot(2, 2, 4);
    %     hold on;    plot(rocket1.t_values, rocket1.x_vel, 'o');
    %     x_vel=own_solver.solve_interpolate(rocket1.t_values,rocket1.x_vel,tt);
    %     % x_vel=own_solver.solve_spline(rocket1.t_values,rocket1.x_vel,tt);
    %     hold on;    plot(tt, x_vel, ':');
    %     x_vel=spline(rocket1.t_values,rocket1.x_vel,tt);
    %     hold on;    plot(tt, x_vel, '-');
%}



f2=figure;

%Solve small rocket best launch time.
little_rockets= createArray(1,number_of_small_rockets,"Rocket");
land_points=zeros(2,number_of_small_rockets);
t_interval=[0,rocket1.t_values(rocket1_max_height.index+1)];
prev_time_to_fire = 0;

disp(newline)
disp("Calculating when to fire small rocket...")
for iteration =1:5
    %Interpolate start conditions for small rocket.
    small_rocket_start_times=linspace(t_interval(1),t_interval(2),number_of_small_rockets);

    little_rocket_y_start=own_solver.solve_interpolate(rocket1.t_values,rocket1.y_pos,small_rocket_start_times);
    little_rocket_x_start=own_solver.solve_interpolate(rocket1.t_values,rocket1.x_pos,small_rocket_start_times);
    little_rocket_y_vel=own_solver.solve_interpolate(rocket1.t_values,rocket1.y_vel,small_rocket_start_times);
    little_rocket_x_vel=own_solver.solve_interpolate(rocket1.t_values,rocket1.x_vel,small_rocket_start_times);
    %{
        % little_rocket_y_start=own_solver.solve_spline(rocket1.t_values,rocket1.y_pos,small_rocket_start_times);
        % little_rocket_x_start=own_solver.solve_spline(rocket1.t_values,rocket1.x_pos,small_rocket_start_times);
        % little_rocket_y_vel=own_solver.solve_spline(rocket1.t_values,rocket1.y_vel,small_rocket_start_times);
        % little_rocket_x_vel=own_solver.solve_spline(rocket1.t_values,rocket1.x_vel,small_rocket_start_times);
        % little_rocket_y_start=spline(rocket1.t_values,rocket1.y_pos,small_rocket_start_times);
        % little_rocket_x_start=spline(rocket1.t_values,rocket1.x_pos,small_rocket_start_times);
        % little_rocket_y_vel=spline(rocket1.t_values,rocket1.y_vel,small_rocket_start_times);
        % little_rocket_x_vel=spline(rocket1.t_values,rocket1.x_vel,small_rocket_start_times);
    %}

    %Solve and plot small rocket trajectories
    for i = 1:number_of_small_rockets
        start_x_pos=little_rocket_x_start(i);
        start_y_pos=little_rocket_y_start(i);
        start_x_vel=little_rocket_x_vel(i);
        start_y_vel=little_rocket_y_vel(i);

        %Solve trajectory and get landing points of small rockets
        little_rockets(i)=Rocket(start_x_pos,start_y_pos,start_x_vel,start_y_vel,rocket2_fuel,rocket2_mass,rocket2_burn_time,rocket2_air_resistance,rocket2_burn_force).solve_trajectory(rocket2_trajectory_end_time, tolerance);
        land_points(:,i)=little_rockets(i).get_land_point().point;
    end
    
    %Find smaller interval to launch rockets
    best_times=find_three_best_candidates(small_rocket_start_times,land_points(1,:));
    t_interval(1)=min(best_times.t_values);
    t_interval(2)=max(best_times.t_values);

    %Finsh by fitting a second degree polynom and solving for the optimal launch time
    max_x_distance_small_rocket_polynom = own_solver.solve_polyfit(best_times.t_values,land_points(1,best_times.indexes),2);
    p = max_x_distance_small_rocket_polynom;
    time_to_fire = -p(2)/(2*p(1));

    E_trunk = abs(time_to_fire - prev_time_to_fire);
    disp("Time to fire: " + sprintf('%.10f', time_to_fire) + "   trajectory E_trunk: " + sprintf('%.10f', little_rockets(best_times.indexes(1)).get_trajectory_E_trunk()) + "   E_trunk: " + sprintf('%.10f', E_trunk) + "   global E_trunk: " + sprintf('%.10f', little_rockets(best_times.indexes(1)).get_trajectory_E_trunk() + E_trunk) + "      (Distance reached: " + sprintf('%.10f', own_solver.solve_polyval(p, time_to_fire)) + ")");

    prev_time_to_fire = time_to_fire;

    %Plot landing points of small rockets against launch time
    hold on;    plot(small_rocket_start_times,land_points(1,:),'-o')
    %Plot fitted polynom
    tt=linspace(3, 4, 1000);
    hold on;    plot(tt, own_solver.solve_polyval(max_x_distance_small_rocket_polynom, tt), "-.");
    hold on;    plot([time_to_fire, time_to_fire], [-200, 300], "-")
end



%Make plots pretty
title("Distance reached against small rocket start time");
xlabel("Time");
ylabel("Distance reached");
ylim([-20, 400]);






%Plot trajectory for optimal small rocket launch time
f3=figure;
pbaspect([1 1 1]);
max_x_distance_index = best_times.indexes(1);                                                                                      
hold on;    plot(little_rockets(max_x_distance_index).x_pos,little_rockets(max_x_distance_index).y_pos,'*-');       
hold on;    plot([0,300], [0,0]);                                                                                   
hold on;    plot(land_points(1,max_x_distance_index), land_points(2,max_x_distance_index), "o",'MarkerSize',10);     
hold on;    plot(rocket1.x_pos, rocket1.y_pos, ":")


%Make plots pretty
title("Small rocket trajectory");
xlabel("x position");
ylabel("y position");
legend(["Small rocket trajectory", "Ground", "Land point.   x: " + land_points(1,max_x_distance_index) + ",   y:" + land_points(2,max_x_distance_index), "Large rocket trajectory"], 'Location', 'southwest');
