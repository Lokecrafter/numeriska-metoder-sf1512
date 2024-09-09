close all, clear all, clc;

global t_days;
global sun_time_minutes;
t_days = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 366]';
sun_time_minutes = [377, 485, 632, 794, 957, 1083, 1105, 997, 847, 683, 526, 395, 374]';
fine_t_days = linspace(1, 366, 100);

function ret = plot_discrete_data(title_msg)
    global t_days;
    global sun_time_minutes;
    plot(t_days, sun_time_minutes, "o");
    xlabel("Days since new year.");
    ylabel("Sun time in minutes");
    title(title_msg);
    
    hold on;
end

%Linear
subplot(3, 3, 1);
plot_discrete_data("Linear interpolation");
plot(t_days, sun_time_minutes, "-");

%One polynomial
subplot(3, 3, 2);

A = [t_days .^ 0, t_days .^ 1, t_days .^ 2, t_days .^ 3, t_days .^ 4, t_days .^ 5, t_days .^ 6, t_days .^ 7, t_days .^ 8, t_days .^ 9, t_days .^ 10, t_days .^ 11, t_days .^ 12];
c = A\sun_time_minutes;

polynomial_interpolation = polyval(c', fine_t_days);

plot_discrete_data("Polynomial interpolation");
plot(fine_t_days, polynomial_interpolation);

%Splines interpolation
subplot(3, 3, 8);
sun_time_spline = spline(t_days, sun_time_minutes, fine_t_days);
plot_discrete_data("Linear interpolation");
plot(fine_t_days, sun_time_spline, "-");
