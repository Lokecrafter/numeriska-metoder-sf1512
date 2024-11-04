clear all; clc; clf;



% classdef BasicClass
%     properties
%        Value {mustBeNumeric}
%     end
%     methods
%        function r = roundOff(obj)
%           r = round([obj.Value],2);
%        end
%        function r = multiplyBy(obj,n)
%           r = [obj.Value]*n;
%        end
%     end
% end


classdef Rocket
    properties
       x_pos
       y_pos
       x_vel
       y_vel
       mass
       min_mass
       mass_velocity
       air_resistance
    end
    methods
       function r = roundOff(obj)
          r = round([obj.Value],2);
       end
       function r = multiplyBy(obj,n)
          r = [obj.Value]*n;
       end
    end
end