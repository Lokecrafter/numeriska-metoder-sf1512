close all, clear all, clc;

global data_t_days;
global data_sun_time_minutes;
data_t_days = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 366]';
data_sun_time_minutes = [377, 485, 632, 794, 957, 1083, 1105, 997, 847, 683, 526, 395, 374]';
fine_t_days = linspace(min(data_t_days), max(data_t_days), 100);

saved_interpolations_matrix = [];

function ret = plot_discrete_data(title_msg)
    global data_t_days;
    global data_sun_time_minutes;
    plot(data_t_days, data_sun_time_minutes, "o");
    xlabel("Days since new year");
    ylabel("Sun time in minutes");
    title(title_msg);
    xlim([min(data_t_days), max(data_t_days)]);
    ylim([0, 1200]);
    
    hold on;
end

%12-deg polynomial
subplot(3, 3, 1);
[p,~,centering] = polyfit(data_t_days, data_sun_time_minutes, length(data_t_days) - 1);
interpolation = polyval(p, fine_t_days,[],centering);
saved_interpolations_matrix = [saved_interpolations_matrix ; interpolation];
plot_discrete_data("12-deg polynomial interpolation");
plot(fine_t_days, interpolation);


%Linear
subplot(3, 3, 2);
interpolation = interp1(data_t_days, data_sun_time_minutes, fine_t_days);
saved_interpolations_matrix = [saved_interpolations_matrix ; interpolation];
plot_discrete_data("Linear interpolation");
plot(fine_t_days, interpolation, "-");


%Splines interpolation
subplot(3, 3, 3);
interpolation = spline(data_t_days, data_sun_time_minutes, fine_t_days);
saved_interpolations_matrix = [saved_interpolations_matrix ; interpolation];
plot_discrete_data("Splines interpolation");
plot(fine_t_days, interpolation, "-");


%June-August polynomial
subplot(3, 3, 4);
[p,~,centering] = polyfit(data_t_days(6:8), data_sun_time_minutes(6:8), 2);
interpolation = polyval(p, fine_t_days,[],centering);
saved_interpolations_matrix = [saved_interpolations_matrix ; interpolation];
plot_discrete_data("2-deg polynomial interpolation (june-august)");
plot(fine_t_days, interpolation);


%April-September polynomial
subplot(3, 3, 5);
[p,~,centering] = polyfit(data_t_days(4:9), data_sun_time_minutes(4:9), 2);
interpolation = polyval(p, fine_t_days,[],centering);
saved_interpolations_matrix = [saved_interpolations_matrix ; interpolation];
plot_discrete_data("2-deg polynomial interpolation (april-september)");
plot(fine_t_days, interpolation);


%Januari-(31th December) polynomial
subplot(3, 3, 6);
[p,~,centering] = polyfit(data_t_days, data_sun_time_minutes, 2);
interpolation = polyval(p, fine_t_days,[],centering);
saved_interpolations_matrix = [saved_interpolations_matrix ; interpolation];
plot_discrete_data("2-deg polynomial interpolation (Januari-(31th December))");
plot(fine_t_days, interpolation);




%Januari-(31th December) polynomial
subplot(3, 3, 6);
[p,~,centering] = polyfit(data_t_days, data_sun_time_minutes, 2);
interpolation = polyval(p, fine_t_days,[],centering);
saved_interpolations_matrix = [saved_interpolations_matrix ; interpolation];
plot_discrete_data("2-deg polynomial interpolation (Januari-(31th December))");
plot(fine_t_days, interpolation);