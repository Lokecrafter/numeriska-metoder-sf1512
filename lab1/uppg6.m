close all, clear all, clc;

%integrand = @(x) (1 - exp(-(x/3).^3)) ./ (4*x.^3);
integrand = @(x) (1 - exp(-(x/3).^3)) ./ (4*x.^3);



%x = linspace(0, 10^(-4), 10000);

%plot(x, integrand(x));
%hold on;
%plot([0,0], [0,10^(-4)]);


%Integranden oscillerar jättemycket. Detta beror troligen på trunkeringsfel i datorn eftersom värdet är så litet på täljaren och nämnaren i integranden.

%Vi integrerar till x-värdet B. Vi hittar en integrand som är större än eller lika med ursprungliga integranden och sen beräknar vi felet vi får av att klippa av integralen.
%Vi kan inte gå från noll heller.

%upper_limit = 10;
%lower_limit = 10^-(10);

%x = linspace(lower_limit, upper_limit, 100);


%integrand = @(x) 100*x.^5;

% FLERA TRAPETSER
a= 0.001; b= 22361; n= 2;
t=0;
for i = 1:10;
    h = (b - a)/n;
    x = a + h * [0:n];
    y = integrand(x);
    t(i) = h * (sum(y) - (y(1) + y(n+1)) / 2);
    n = 2 * n;
end;
t=t';

t