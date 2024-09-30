%Solving differential equation for pendulum.

%Own solver
% close all; clear all; clc;

% start_y = 5 * pi / 7;
% start_y_prim = 0.9;
% g = 9.82;
% L = 2.35;

% acceleration = @(y) - sin(y) * g / L;

% y_vect = [start_y; start_y_prim];
% disp(y_vect)

% h = 0.1;

% next_y_vect = y_vect + h * [y_vect(2); acceleration(y_vect(1))];
% disp(next_y_vect)

g = 9.82;
L = 2.35;

%fi  = fi  + h fi'
%fi' = fi' + h fi''

odefun = @y [1; ];
tspan = [0, 1];
y0 = [5 * pi / 7; 0.9];

[t,y] = ode45(odefun,tspan,y0);