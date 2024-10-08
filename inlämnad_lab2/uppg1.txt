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

%-----------------------------------------------

% DELUPPGIFT b

%-----------------------------------------------

clear all; clc; clf;

%fi  = fi  + h fi'
%fi' = fi' + h fi''

T = 12;
L = 2.35;
% Adjusted for Deluppgift e
%L = 2.5;
G = 9.82;
h = 0.1;

odefun = @(t,y) [y(2); - G * sin(y(1))/ L ]; % gives the vector "u'_n"
y_bis =@(t,y) - G * sin(y(1)) / L ;
tspan = [0, T];
y0 = [5 * pi / 7; 0.9];

    
% Using the ode45 function to calculate
[t,y] = ode45(odefun,tspan,y0);
% disp([t, y])
% plot(t, y, '-o');
% xlabel = t;
% ylabel = y;
% hold on
% plot(tspan, zeros(size(tspan)), '-');

%-----------------------------------------------

% DELUPPGIFT c

%-----------------------------------------------

% Given function for animation
function anim(tut, fiut, L)
    for i=1:length(tut)-1
        x0 = L * sin(fiut(i)); 
        y0 = - L * cos(fiut(i));
        plot([0,x0],[0,y0],'-o')
        axis('equal')
        axis([-1 1 -1 0]*1.2*L)
        drawnow
        pause(tut(i+1)-tut(i))
    end
end

% anim(t, y(:,1), L) % y(:,1) in the y-matrix, all rows, in column 1

%-----------------------------------------------

% DELUPPGIFT d, e

%-----------------------------------------------

% Över vilket eller vilka intervall bör interpolationen gå? 
% Är långt eller kort intervall bäst? 
% Vilken grad väljer ni? 
% Hur många polynom? 
% Hur skattar ni osäkerheten? 
% (Jämför gärna era beslut här med era slutsatser från Lab1!)

% Splines interpol. 
tt = linspace(tspan(1), tspan(2), 10000);
interpolation = spline(t, y(:,1), tt);
plot(tt, interpolation, "-");
hold on

% Calculates index for a given t-value.
index =@(x) x * length(interpolation) / tspan(2);

% Peak 1
disp(['Peak 1 indexes: ', num2str(index(4)), ' ', num2str(index(6))])
% Rounded manually

% Peak 2
disp(['Peak 1 indexes: ', num2str(index(9)), ' ', num2str(index(11))])
% Rounded manually

% Manually input for interval indexes
[max1, index1] = max(interpolation(3333:5000));
[max2, index2] = max(interpolation(7500:9166));
disp(['Max 1 and max 2: ', num2str(max1), ' ', num2str(max2)]) % check to have max1 and max2 be ish equal
disp(['Index for max 1 and max 2: ', num2str(index1), ' ', num2str(index2)]) % Manual check for index

tt1 = tt(3333 + index1); % Max1's t-value based on manual inputs for indexes
tt2 = tt(7500 + index2); % Max2's -II-
disp(['Max1, and 2s t-value based on manual inputs for indexes: ' num2str(tt1), ' ', num2str(tt2)]) % con
plot([tt1, tt2], [max1, max2], 'o'); % marks our max-values
period_time = tt2 - tt1; 
disp(['Period time: ', num2str(period_time), ' seconds'])