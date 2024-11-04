classdef Rocket
    properties
       x_pos            {mustBeNumeric}
       y_pos            {mustBeNumeric}
       x_vel            {mustBeNumeric}
       y_vel            {mustBeNumeric}
       fuel_mass        {mustBeNumeric}
       body_mass        {mustBeNumeric}
       burn_time        {mustBeNumeric}
       air_resistance   {mustBeNumeric}
       force            {mustBeNumeric}
    end
    methods
        function obj=Rocket(x_pos0,y_pos0,x_vel0,y_vel0,fuel_mass_g, body_mass_g,new_burn_time_s,new_air_resistance,new_force)
            if nargin == 9
                obj.x_pos = x_pos0;
                obj.y_pos = y_pos0;
                obj.x_vel = x_vel0;
                obj.y_vel = y_vel0;
                obj.fuel_mass = fuel_mass_g;
                obj.body_mass = body_mass_g;
                obj.burn_time = new_burn_time_s;
                obj.air_resistance = new_air_resistance;
                obj.force=new_force;
            end
        end
    end
end