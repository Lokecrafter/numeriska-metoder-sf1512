% lab 2 upg 3

%-----------------------------------------------

% Just det här differentialekvationsproblemet i uppgift 2 är väldigt enkelt att lösa med
% inskjutningsmetoden. Skriv ett Matlab-program som skattar temperaturen i x = 1.40
% med inskjutningsmetoden och Matlabs ode45 och fzero. Glöm inte felskattningen.

%-----------------------------------------------

clear all; close all; clc; clf;

L = 3.40; % [m]
T0 = 300; % [K]
TL = 450; % [K] 

% Söker x = 1.40 [m]

n = 34; % Nr of steps (adjusted to "hit" x = 1.40)

tolerance = 1e-4;

prev_prev_T_index = 2;
prev_T_index = 1;



Q  =@(x) 285 .* exp(- ((x - L .* 0.5).^2));
    
C1 =@(x) (3 + x .* 1 ./ 6) .* (1 ./ (h.^2)) .* (-2);                    % For T_i-1
C2 =@(x) (3 + x .* 1 ./ 6) .* (1 ./ (h.^2)) - (1 ./ (6 .* 2 .* h));     % For T_i
C3 =@(x) (1 ./ (6 .* 2 .* h) + (3 + x .* 1 ./ 6) .* (1 ./ (h.^2)));     % For T_i+1

disp('For x = 1.40, T(x) = : ')
disp('   T          prev_T prev_prev_T E_trunk convergence')