% Uppgift 8 - Numerisk integration: Eulers metod



%-----------------------------------------------
%-----------------------------------------------
%-----------------------------------------------
% INGEN KOD BEHÖVER LADDAS UPP FÖR DENNA UPPGIFT
%-----------------------------------------------
%-----------------------------------------------
%-----------------------------------------------



clear; clf; clc;

% Eulers metod är en grundläggande teknik för numerisk integration, 
% specifikt för att lösa ordinära differentialekvationer (ODE).

% Metoden används för att lösa ekv. av typen dy/dt = f(t,y) med
% initialvillkoret y(t_0) = y_0. 
% I denna uppgift är det dy/dx=f(y,x) och vi har villkoret y(x_0)=y_0.

% Steg:
% 1.   Antag ovanstående villkor för ekvationen.
% 2.   Välj steglängd h.
% 3.   Iterera med Eulers (iterations)metod.

%-----------------------------------------------

% DELUPPGIFT A

%-----------------------------------------------

% 1.   dy/dx och initialvillkoret:
%---------------------------------
f=@(fx, fy) -(1/6 + (pi.*sin(pi.*fx))/(1.6 - cos(pi.*fx))).*fy;
x_0 = 0; % x initialvillkorsvärde.
y_0 = 2.5; % y(x) initialvillkorsvärde, dvs y_0 = y(0) i detta fall.

% 2.   Väljer steglängd h:
%-------------------------
h = 0.5;

% 3. Eulers metod (startkrav):
%-----------------------------
x_max = 6; % Givet i upg. Används som maxvärdet för intervallet nedan.
x = x_0:h:x_max; % Intervallet för x vi vill jobba på.
y = zeros(size(x)); % Skapar y-vektor lika stor som x-vektorn som 
% skapades i raden ovanför. zeros gör att hela vektorn fylls av
% nollor (0) vilket även säkerställer att ingen tidigare data kan råka
% ligga i vår vektor.
y(1)=y_0; % Vi ger första värdet för y-vektorn 
% som vi vet är initialvillkorets värde givet tidigare.


% 3. Eulers metod (iterationsloop):
prev_y_6 = 2.5;
for i = 1:15
    x = x_0:h:x_max; % Intervallet för x vi vill jobba på.
    y = zeros(size(x));
    y(1) = y_0;
    %----------------------------------
    for a = 1:length(x) - 1 % Längden av x fast vi vill inte göra en 
        % gång för många så det blir length(x)-1. n är vår loop-variabel.
        y(a+1) = y(a) + h * f(x(a), y(a)); % Beräknar nästa y-värde med 
        % Eulers metod! So simple! <3 :D
    end
    
    %disp("h = " + h + "   y(6) = " + y(end) + "    E_trunk = " + abs(y(end) - prev_y_6) + "   y(6) + E_trunk = " + (y(end) + abs(y(end) - prev_y_6)))
    disp("h = " + h + "   y(6) = " + y(end) + "   y(6) + E_trunk = " + (y(end) + abs(y(end) - prev_y_6)) + "   E_trunk = " + (abs(y(end) - prev_y_6)))
    %disp("E_trunk = " + (abs(y(end) - prev_y_6)))
    %disp((y(end) + abs(y(end) - prev_y_6)));
    prev_y_6 = y(end);
    h = h *0.5;
    
end
plot(x,y,'-o'); % '-o' : 
% '-' Skapar en solid linje som förbinder datapunkterna. 
% 'o': Lägger till cirkelmarkörer vid varje datapunkt.
xlabel('x');
ylabel('y');
grid on;
title('Eulers metod för Upg 8')


hold on;
%ode45(f, [0,6], 2.5)