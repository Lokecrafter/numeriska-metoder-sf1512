clear all; clc; close all;
import Rocket.*;

global rocket1;
global rocket2;

rocket1=Rocket(0,0,0,0,45,20,3,0.01,10);
rocket2=Rocket(0,0,0,0,12,5,4,0.001,3);

u_0=[rocket1.x_pos;rocket1.y_pos;rocket1.x_vel;rocket1.y_vel];

function ret=odefun(t,u)
    global rocket1;
    global rocket2;

    ret=[0;0;0;0];
    ret(1)=u(3);
    ret(2)=u(4);

    V=hypot(u(3),u(4));
    K=rocket1.air_resistance;
    F=rocket1.force*(t<=rocket1.burn_time); %istället för en if-sats
    mass=rocket1.body_mass+rocket1.fuel_mass*max(0, (rocket1.burn_time-t)/rocket1.burn_time);
    mass=mass/1000; %massan i kilo
    angle=0;
    if u(3)==0 && u(4) == 0 
        angle=deg2rad(80);
    else
        angle=atan2(u(4),u(3));
    end

    ret(3)=F*cos(angle)-K*u(3)*V;
    ret(4)=F*sin(angle)-K*u(4)*V-mass*9.82;
    ret(3)=ret(3)/mass;
    ret(4)=ret(4)/mass;
end

tspan=[0,20];

result=ode45(@odefun,tspan,u_0);
%plot(result.x,result.y,'o-')
plot(result.y(1,:),result.y(2,:),'-o') %y(1) är raketens x-position, y(2) är raketens y-position
hold on 
plot([0,70],[0,0])