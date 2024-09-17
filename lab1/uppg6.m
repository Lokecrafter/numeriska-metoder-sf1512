close all, clear all, clc;

integrand = @(x) (1 - exp(-(x/3).^3)) ./ (4*x.^3);



x = linspace(0, 10^(-4), 10000);

plot(x, integrand(x));