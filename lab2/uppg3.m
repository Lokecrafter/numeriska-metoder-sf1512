% lab 2 upg 3

%-----------------------------------------------

% Just det här differentialekvationsproblemet i uppgift 2 är väldigt enkelt att lösa med
% inskjutningsmetoden. Skriv ett Matlab-program som skattar temperaturen i x = 1.40
% med inskjutningsmetoden och Matlabs ode45 och fzero. Glöm inte felskattningen.

%-----------------------------------------------

clear all; close all; clc;
global L;
global T0;
global TL;
L = 3.40; % [m]
T0 = 300; % [K]
TL = 450; % [K] 

% Söker x = 1.40 [m]

% n = 34; % Nr of steps (adjusted to "hit" x = 1.40)

tolerance = 1e-4;

% prev_prev_T_index = 2;
% prev_T_index = 1;


% C1 =@(x) (3 + x .* 1 ./ 6) .* (1 ./ (h.^2)) .* (-2);                    % For T_i-1
% C2 =@(x) (3 + x .* 1 ./ 6) .* (1 ./ (h.^2)) - (1 ./ (6 .* 2 .* h));     % For T_i
% C3 =@(x) (1 ./ (6 .* 2 .* h) + (3 + x .* 1 ./ 6) .* (1 ./ (h.^2)));     % For T_i+1

% disp('For x = 1.40, T(x) = : ')
% disp('   T          prev_T prev_prev_T E_trunk convergence')

%-----------------------------------------------
%
% Gives T_prim_start_guess
%
%-----------------------------------------------

function error = error_function(T_prim_guess)
    global L;
    global T0;
    global TL;

    Q  =@(x) 285 .* exp(- ((x - L .* 0.5).^2));
    
    odefun =@(t, y) [y(2); ( - Q(t) - (y(2) / 6 )) / (3 + t / 6)]; 
    x_span = [0,L];
    y0 = [T0;T_prim_guess];
    
    result = ode45(odefun, x_span, y0);
    
    error = result.y(1,end) - TL; % Error - checking if the guessed T_prim_guess gives the right value for a sequence which gives TL. 
end

xx = linspace(0, 1000, 1000);
yy = zeros(size(xx));

for i = 1:length(yy)
    yy(i) = error_function(xx(i));
end

f1 = figure;
plot(xx, yy, '-')
hold on
plot([0,1000],[0,0])
hold on

T_prim_start_guess = 128; % For fzero

T_prim_start_value = fzero(@error_function, T_prim_start_guess); % Matlab's sekantmetod
disp("Solved T_prim-value with Sekantmetoden: " + T_prim_start_value + " K/m")
plot(T_prim_start_value, error_function(T_prim_start_value), 'o') % Visual control

%-----------------------------------------------
%
% T(x)
%
%-----------------------------------------------

Q  =@(x) 285 .* exp(- ((x - L .* 0.5).^2));

odefun =@(t, y) [y(2); ( - Q(t) - (y(2) / 6 )) / (3 + t / 6)]; 
x_span = [0,L];
y0 = [T0;T_prim_start_value];

result = ode45(odefun, x_span, y0);

f2 = figure;
plot(result.x, result.y, '-o');

T_interpol_for_x_value = spline(result.x, result.y, 1.4); % ans in T(1.4) and T'(1.4) (red)
disp("T(1.4) = " + T_interpol_for_x_value(1) + " K, T'(1.4) = " + T_interpol_for_x_value(2) + " K/m")
