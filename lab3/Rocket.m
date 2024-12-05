
classdef Rocket
    properties
        t_values 
        x_pos            
        y_pos            
        x_vel            
        y_vel            
        fuel_mass        {mustBeNumeric}
        body_mass        {mustBeNumeric}
        burn_time        {mustBeNumeric}
        air_resistance   {mustBeNumeric}
        force            {mustBeNumeric}
        solver
        E_trunk_trajectory {mustBeNumeric}
    end
    methods
        function obj=Rocket(x_pos0,y_pos0,x_vel0,y_vel0,fuel_mass_g, body_mass_g,new_burn_time_s,new_air_resistance,new_force)
            import Solvers.*;
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
                obj.solver = Solvers(true);
            end
        end
        function ret=solve_trajectory(obj,end_time,tolerance)
            function ret=odefun(t,u)
                ret=[0;0;0;0];
                ret(1)=u(3);
                ret(2)=u(4);
                
                V=hypot(u(3),u(4));
                K=obj.air_resistance;
                F=obj.force*(t<obj.burn_time); %instead of an if-sats
                mass=obj.body_mass+obj.fuel_mass*max(0, (obj.burn_time-t)/obj.burn_time);
                mass=mass/1000; %mass in kilo
                angle=0;
                if u(3)==0 && u(4) == 0 
                    angle=deg2rad(80);
                else
                    angle=atan2(u(4),u(3));
                end
                
                ret(3)=F*cos(angle)-K*u(3)*V;
                ret(4)=F*sin(angle)-K*u(4)*V-mass*9.82;
                ret(3)=ret(3)/mass;
                ret(4)=ret(4)/mass;
            end
                 
            obj.E_trunk_trajectory = 0;

            %solve interval for rocket burn time
            u_0=[obj.x_pos(end);obj.y_pos(end);obj.x_vel(end);obj.y_vel(end)];
            tspan=[0,obj.burn_time];
            result=obj.solver.solve_ode45(@odefun,tspan,u_0,50,tolerance/3);
            
            obj.t_values=result.x;
            obj.x_pos=result.y(1,:);
            obj.y_pos=result.y(2,:);
            obj.x_vel=result.y(3,:);
            obj.y_vel=result.y(4,:);
            obj.E_trunk_trajectory = obj.E_trunk_trajectory + result.E_trunk;
            
            
            %solve interval for rocket slowing time
            u_0=[obj.x_pos(end);obj.y_pos(end);obj.x_vel(end);obj.y_vel(end)];
            tspan=[obj.burn_time,obj.burn_time+2];
            result=obj.solver.solve_ode45(@odefun,tspan,u_0,10,tolerance/3);
            
            obj.t_values=[obj.t_values,result.x(2:end)];
            obj.x_pos=[obj.x_pos,result.y(1,2:end)];
            obj.y_pos=[obj.y_pos,result.y(2,2:end)];
            obj.x_vel=[obj.x_vel,result.y(3,2:end)];
            obj.y_vel=[obj.y_vel,result.y(4,2:end)];
            obj.E_trunk_trajectory = obj.E_trunk_trajectory + result.E_trunk;

            
            
            %solve interval for rocket falling time
            u_0=[obj.x_pos(end);obj.y_pos(end);obj.x_vel(end);obj.y_vel(end)];
            tspan=[obj.burn_time+2,end_time];
            result=obj.solver.solve_ode45(@odefun,tspan,u_0,10,tolerance/3);
            
            obj.t_values=[obj.t_values,result.x(2:end)];
            obj.x_pos=[obj.x_pos,result.y(1,2:end)];
            obj.y_pos=[obj.y_pos,result.y(2,2:end)];
            obj.x_vel=[obj.x_vel,result.y(3,2:end)];
            obj.y_vel=[obj.y_vel,result.y(4,2:end)];
            obj.E_trunk_trajectory = obj.E_trunk_trajectory + result.E_trunk;
            
            ret=obj;
        end
        function ret=get_land_point(obj)
            prev_index=10;
            for i = prev_index:length(obj.y_pos)
                prev_sign=sign(obj.y_pos(prev_index));
                current_sign=sign(obj.y_pos(i));
                if current_sign ~= prev_sign
                    break
                end
                prev_index=i;
            end
            %calculates when the rocket crosses the x-axis
            pt1=[obj.x_pos(prev_index-1);obj.y_pos(prev_index-1)];
            pt2=[obj.x_pos(prev_index);obj.y_pos(prev_index)];
            direction=pt2-pt1;
            
            land_point=pt1-pt1(2).*[direction(1)/direction(2);1]; 
            E_trunk = max(abs(land_point(1) - pt1(1)), abs(land_point(1) - pt2(1)));

            ret.point=land_point;
            ret.E_trunk = E_trunk;
            ret.glob_E_trunk = E_trunk + obj.E_trunk_trajectory;
        end
        function ret=get_highest_point(obj)
            [~,max_height_index]=max(obj.y_pos);

            xx = [obj.x_pos(max_height_index-1:max_height_index+1)];
            yy = [obj.y_pos(max_height_index-1:max_height_index+1)];

            p = obj.solver.solve_polyfit(xx,yy,2);


            x_max_height = -p(2)/(2*p(1));
            y_max_height = obj.solver.solve_polyval(p,x_max_height);

            x = x_max_height;
            y = y_max_height;
            index = max_height_index;
            E_trunk = norm([x_max_height; y_max_height] - [obj.x_pos(max_height_index); obj.y_pos(max_height_index)]);

            ret.point = [x,y];
            ret.index = index;
            ret.E_trunk = E_trunk;
            ret.glob_E_trunk = E_trunk + obj.E_trunk_trajectory;
        end
        function ret=get_trajectory_E_trunk(obj)
            ret = obj.E_trunk_trajectory;
        end
    end
end