classdef Rocket
    properties
       x_pos            {mustBeNumeric}
       y_pos            {mustBeNumeric}
       x_vel            {mustBeNumeric}
       y_vel            {mustBeNumeric}
       mass             {mustBeNumeric}
       min_mass         {mustBeNumeric}
       mass_velocity    {mustBeNumeric}
       air_resistance   {mustBeNumeric}
    end
    methods
        function obj=Rocket(x_pos0,y_pos0,x_vel0,y_vel0,mass0_g, new_min_mass_g,new_mass_velocity_g_per_s,new_air_resistance)
            if nargin == 8
                obj.x_pos = x_pos0;
                obj.y_pos = y_pos0;
                obj.x_vel = x_vel0;
                obj.y_vel = y_vel0;
                obj.mass = mass0_g;
                obj.min_mass = new_min_mass_g;
                obj.mass_velocity = new_mass_velocity_g_per_s;
                obj.air_resistance = new_air_resistance;
            end
        end
    end
end