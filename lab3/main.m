clear all; clc; close all;
import Rocket.*;
f1=figure;


rocket1=Rocket(0,0,0,0,45,20,3,0.01,10);

%u_0=[rocket1.x_pos;rocket1.y_pos;rocket1.x_vel;rocket1.y_vel];

%function ret=odefun(t,u)
%    global rocket1;
%    global rocket2;
%
%    ret=[0;0;0;0];
%    ret(1)=u(3);
%    ret(2)=u(4);
%
%    V=hypot(u(3),u(4));
%    K=rocket1.air_resistance;
%    F=rocket1.force*(t<=rocket1.burn_time); %istället för en if-sats
%    mass=rocket1.body_mass+rocket1.fuel_mass*max(0, (rocket1.burn_time-t)/rocket1.burn_time);
%    mass=mass/1000; %massan i kilo
%    angle=0;
%    if u(3)==0 && u(4) == 0 
%        angle=deg2rad(80);
%    else
%        angle=atan2(u(4),u(3));
%    end
%
%    ret(3)=F*cos(angle)-K*u(3)*V;
%    ret(4)=F*sin(angle)-K*u(4)*V-mass*9.82;
%    ret(3)=ret(3)/mass;
%    ret(4)=ret(4)/mass;
%end

%tspan=[0,20];

%result=ode45(@odefun,tspan,u_0);
%plot(result.x,result.y,'o-')
rocket1=rocket1.solve_trajectory(20);
%plot(result.y(1,:),result.y(2,:),'-o') %y(1) är raketens x-position, y(2) är raketens y-position
plot(rocket1.x_pos,rocket1.y_pos,'-o') %y(1) är raketens x-position, y(2) är raketens y-position
hold on 
plot([0,70],[0,0])

%rocket1.t_values=result.x;
%rocket1.x_pos=result.y(1,:);
%rocket1.y_pos=result.y(2,:);
%rocket1.x_vel=result.y(3,:);
%rocket1.y_vel=result.y(4,:);

[~,max_height_index]=max(rocket1.y_pos);

xx=[rocket1.x_pos(max_height_index-1:max_height_index+1)];
yy=[rocket1.y_pos(max_height_index-1:max_height_index+1)];

%[p,~,mu]=polyfit(xx,yy,2);
p=polyfit(xx,yy,2);
xxx=linspace(60,70,100);

hold on

%plot(xxx,polyval(p,xxx,[],mu))
hold on
plot(xxx,polyval(p,xxx))

x_max_height=-p(2)/(2*p(1))
y_max_height=polyval(p,x_max_height) %SVAR på uppgift a
hold on
plot(x_max_height,y_max_height,'o')

prev_index=10;
for i = prev_index:length(rocket1.y_pos)
    prev_sign=sign(rocket1.y_pos(prev_index));
    current_sign=sign(rocket1.y_pos(i));
    if current_sign ~= prev_sign
        break
    end
    prev_index=i;
end


%beräknar skärningspunkt med x-axeln
pt1=[rocket1.x_pos(prev_index);rocket1.y_pos(prev_index)];
pt2=[rocket1.x_pos(prev_index+1);rocket1.y_pos(prev_index+1)];
direction=pt2-pt1;

land_point=pt1-pt1(2).*[direction(1)/direction(2);1]; %SVAR på uppgift b

hold on
plot(land_point(1),land_point(2),'o')
disp(land_point)

n=4;
little_rockets= createArray(1,n,"Rocket");

tt=linspace(0,rocket1.t_values(max_height_index+1),n);
little_rocket_y_start=spline(rocket1.t_values,rocket1.y_pos,tt);
little_rocket_x_start=spline(rocket1.t_values,rocket1.x_pos,tt);
little_rocket_y_vel=spline(rocket1.t_values,rocket1.y_vel,tt);
little_rocket_x_vel=spline(rocket1.t_values,rocket1.x_vel,tt);

f2=figure;
for i = length(tt)
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
little_rockets(1).x_pos
little_rockets(2).x_pos
little_rockets(3).x_pos
little_rockets(4).x_pos